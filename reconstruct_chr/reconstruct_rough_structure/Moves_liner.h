#ifndef INVERSE_3D_NEW_MOVES_H
#define INVERSE_3D_NEW_MOVES_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <mutex>
#include <chrono>
#include "Energy_changes_liner.h"

using namespace std;
using namespace Eigen;

int reduce_count=0,add_count=0,kink_count=0,loop_count=0,skip_count=0;

uniform_int_distribution<int> unidir(0,2);
uniform_int_distribution<int> unidir_loop(1,5);
// uniform_int_distribution<int> uniIndex(0,tads_bound.size()-1);

mutex tad_change_mutex;//lock when change tad_change_flags
std::condition_variable cv;
vector<condition_variable> cvs(48);
unordered_map<int,int> tad_change_flags={{0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,0},{7,0},{8,0},{9,0},\
                {10,0},{11,0},{12,0},{13,0},{14,0},{15,0},{16,0},{17,0},{18,0},{19,0},\
                {20,0},{21,0},{22,0},{23,0},{24,0},{25,0},{26,0},{27,0},{28,0},{29,0},\
                {30,0},{31,0},{32,0},{33,0},{34,0},{35,0},{36,0},{37,0},{38,0},{39,0},\
                {40,0},{41,0},{42,0},{43,0},{44,0},{45,0},{46,0},{47,0}};

//---------------------------accept move---------------------------------

bool accept_move(float delta_E){//如果能量变小了，接受Move。如果能量变大了，以一定概率接受Move
    return ((delta_E <=0)||(unif(gen) < exp(-delta_E)));//return 1接受 0拒绝
}

//---------------------------boundary---------------------------------

bool check_boundary(Eigen::Vector3i prop_move1) {
    //find if point in sphere
    if (pow(prop_move1[0],2)+pow(prop_move1[1],2)+pow(prop_move1[2],2)<=pow(radius,2)) { 
        return 1;
    } 
    else {return 0;}
};

// //-------------------apply to the condition which the two site is far----------------------
// void deal_neighbor(vector<Eigen::Vector3i> &polymer ,vector<int> &site_index2position ,vector<vector<int>> &tad_neighbor_contact_new,int number_of_threads,int site,int pol_length,int m,Eigen::Vector3i prop_move){
    
//     //initial tad_neighbor_contact_new
//     if (site%6!=0){return;}

//     int x=prop_move(0);int y=prop_move(1);int z=prop_move(2);
//     find_neighbor(x,y,z,site,m,number_of_threads,polymer,site_index2position,tad_neighbor_contact_new);
//     // cout<<"deal"<<site<<'\t';
//     // for (auto elem:tad_neighbor_contact_new){
//     //     cout<<elem[0]<<'\t';
//     // }
//     // cout<<'\n';
//     return ;

// }

//========================= tad move =======================
float sphere_limit(MatrixXf tads_center,float alpha){
    float sphere_out=0;
    for (int i=0;i<tads_center.rows();i++){
        MatrixXf center=tads_center.row(i);
        float distance2CircleCenter=sqrt(pow(center(0,0),2)+pow(center(0,1),2)+pow(center(0,2),2));
        if (distance2CircleCenter>radius){
            sphere_out+=(distance2CircleCenter-radius)*alpha;
        }
    }
    return sphere_out;
}


void tad_translocation(int thread_num,Eigen::MatrixXi &connected_structure,int &is_move_active,int m){
    vector<Eigen::Matrix<int,1,3>> tad_translate={{2,0,0},{0,2,0},{0,0,2},{-2,0,0},{0,-2,0},{0,0,-2},\
    {1,1,0},{1,-1,0},{1,0,1},{1,0,-1},{-1,1,0},{-1,-1,0},{-1,0,1},{-1,0,-1},\
    {0,1,1},{0,1,-1},{0,-1,1},{0,-1,-1}};

    uniform_int_distribution<int> uniIndex(1,tads_bound.size()-1);
    int tad_index = uniIndex(gen);
    int tad_len=tads_bound.size();
    int first=tads_bound[0][0];
    int start=tads_bound[tad_index][0]-first;
    int end=tads_bound[tad_len-1][1]-first;
    int move_model=uniTranslate(gen);   
    Eigen::MatrixXf tads_center_change(tad_len-tad_index,3);  

    //lock
    mutexs[thread_num].lock();
    Eigen::MatrixXi tadi=connected_structure.block(start,0,end-start+1,3);  
    Eigen::MatrixXi margin=tad_translate[move_model].replicate(tadi.rows(), 1);  
    Eigen::MatrixXi start_bin=connected_structure.row(start-1);

    if(((tadi.row(0)+tad_translate[move_model]-start_bin).cwiseAbs()).sum()<=6){

        tadi+=margin;     

        // get new tad center      
        update_tads_center(tad_index,tadi,tads_center_change);

    }
    else{
        mutexs[thread_num].unlock();
        return;
    }          

    //update contact and location and tad_center
    float sphere_limition=sphere_limit(tads_center_change,alpha)-sphere_limit(tads_center[thread_num].bottomRows(tad_len-tad_index),alpha);
    float energy=delta_E_tads(connected_structure.block(start,0,end-start+1,3),tadi ,tads_center_change,tad_index,thread_num);

    if(accept_move(energy + sphere_limition)){

        //update tad_center
        tads_center[thread_num].bottomRows(tad_len-tad_index)=tads_center_change;  

        update_tads_location(connected_structure.block(start,0,end-start+1,3),tadi ,tad_index,thread_num,m);
        connected_structure.block(start,0,end-start+1,3)=tadi;

        //update tad_data
        for(int index=tad_index;index<tads_bound.size();index++){
            update_tads_center_data(thread_num,index,m);
        }    

    }

    //update m
    threads_m[thread_num]=m;

    //unlock
    mutexs[thread_num].unlock();
}


void tad_rotate(int thread_num,Eigen::MatrixXi &connected_structure,int &is_move_active,int m){
    vector<Eigen::Matrix<int,1,3>> rotate={{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};       
    int tad_len=tads_bound.size();
    uniform_int_distribution<int> uniIndex(0,tad_len-1);
    int tad_index = uniIndex(gen);  
    int first=tads_bound[0][0];
    int last=tads_bound[tad_len-1][1]-first;      
    int dir=unipath(gen);    
    int start=tads_bound[tad_index][0]-first; 
    Eigen::MatrixXf tads_center_change(tad_len-tad_index,3);  


    //rotate
    (mutexs[thread_num]).lock();

    Eigen::MatrixXi rotate_center= connected_structure.row(start); 
    Eigen::MatrixXi rotate_tad=connected_structure.bottomRows(last-start+1);

    Eigen::MatrixXi rotate_mat=rotate[dir].replicate(rotate_tad.rows(), 1);    
    Eigen::MatrixXi margin=rotate_center.replicate(rotate_tad.rows(), 1);
    rotate_tad=(rotate_tad-margin).cwiseProduct(rotate_mat)+margin;

    // get new tad center         
    update_tads_center(tad_index,rotate_tad,tads_center_change);

    //update contact and location   
    float sphere_limition=sphere_limit(tads_center_change,alpha)-sphere_limit(tads_center[thread_num].bottomRows(tad_len-tad_index),alpha);    
    float energy=delta_E_tads(connected_structure.bottomRows(last-start+1),rotate_tad ,tads_center_change,tad_index,thread_num);    
      
    if(accept_move(energy+sphere_limition)){
        
        //update tad_center
        tads_center[thread_num].bottomRows(tad_len-tad_index)=tads_center_change;     

        update_tads_location(connected_structure.bottomRows(last-start+1),rotate_tad ,tad_index,thread_num,m);
        connected_structure.bottomRows(last-start+1)=rotate_tad;

        //update tad_center_data
        for(int index=tad_index;index<tads_bound.size();index++){
            update_tads_center_data(thread_num,index,m);
        }     
    }

    //update m
    threads_m[thread_num]=m;
    //unlock
    (mutexs[thread_num]).unlock();

}

void tad_change(int thread_num,int m){
    int tad_len=tads_bound.size();
    uniform_int_distribution<int> uniIndex(0,tad_len-1);
    //get tad_index
    int tad_index=uniIndex(gen);

    int first=tads_bound[0][0];
    int start=tads_bound[tad_index][0]-first;
    int end=tads_bound[tad_len-1][1]-first;
    int change_thread=uniThread(gen);//交换的thread
    while(thread_num==change_thread){
        change_thread=uniThread(gen);
    }   

    Eigen::MatrixXf tads_center_change1(tad_len-tad_index,3); 
    Eigen::MatrixXf tads_center_change2(tad_len-tad_index,3); 

    // lock and update change_flags
    tad_change_mutex.lock() ;
    if ((tad_change_flags[thread_num]==1) || (tad_change_flags[change_thread]==1)){
        tad_change_mutex.unlock() ;
        return;
    }
    tad_change_flags[thread_num]= 1 ;
    tad_change_flags[change_thread]= 1 ;
    tad_change_mutex.unlock() ;    
  
    (mutexs[change_thread]).lock();

    int m_thread_change= threads_m[change_thread];
    
    //exchange tad
    Eigen::MatrixXi tad1=connected_structure[thread_num].block(start,0,end-start+1,3);
    Eigen::MatrixXi tad2=connected_structure[change_thread].block(start,0,end-start+1,3); 
    
    Eigen::MatrixXi margin_use=tad2.row(0)-tad1.row(0);
    Eigen::MatrixXi margin=(margin_use).replicate(end-start+1 , 1);        

    tad1+=margin;
    tad2-=margin;

    // new tad center         
    update_tads_center(tad_index,tad2,tads_center_change1);
    update_tads_center(tad_index,tad1,tads_center_change2);

    //update contact and location   
    float sphere_limition1=sphere_limit(tads_center_change1,alpha)-sphere_limit(tads_center[thread_num].bottomRows(tad_len-tad_index),alpha);       
    float sphere_limition2=sphere_limit(tads_center_change2,alpha)-sphere_limit(tads_center[change_thread].bottomRows(tad_len-tad_index),alpha);  

    float energy1=delta_E_tads(connected_structure[thread_num].bottomRows(end-start+1),tad2 ,tads_center_change1,tad_index,thread_num); 
    float energy2=delta_E_tads(connected_structure[change_thread].bottomRows(end-start+1),tad1 ,tads_center_change2,tad_index,change_thread);
    
    // float energy1=delta_E(connected_structure[thread_num].bottomRows(end-start+1),tad2 ,tad_index,thread_num); 
    // float energy2=delta_E(connected_structure[change_thread].bottomRows(end-start+1),tad1 ,tad_index,change_thread);


    if(accept_move(energy1+energy2+sphere_limition1+sphere_limition2)){

        //update tad_center
        tads_center[thread_num].bottomRows(tad_len-tad_index)=tads_center_change1;
        tads_center[change_thread].bottomRows(tad_len-tad_index)=tads_center_change2;

        //update location
        update_tads_change_location(connected_structure[thread_num].bottomRows(end-start+1),tad2 ,tad_index,thread_num,m);
        update_tads_change_location(connected_structure[change_thread].bottomRows(end-start+1),tad1 ,tad_index,change_thread,m_thread_change);

        connected_structure[thread_num].block(start,0,end-start+1,3)=tad2;
        connected_structure[change_thread].block(start,0,end-start+1,3)=tad1;

        //update_tads_center data
        for(int index=tad_index;index<tads_bound.size();index++){
            update_tads_center_data(thread_num,index,m);
            update_tads_center_data(change_thread,index,m);
        }            
    }

    //update m
    threads_m[thread_num]=m;    

    //unlock
    (mutexs[change_thread]).unlock();   

    // update tad_match
    for(int i=tad_index;i<tad_len;++i){
        int tmp=tad_match(thread_num,i);
        tad_match(thread_num,i)=tad_match(change_thread,i);
        tad_match(change_thread,i)=tmp;
    }
    
    // update change_flags
    tad_change_mutex.lock();
    tad_change_flags[thread_num]= 0 ;
    tad_change_flags[change_thread]= 0 ;
    tad_change_mutex.unlock();           
}



#endif //INVERSE_3D_NEW_MOVES_H