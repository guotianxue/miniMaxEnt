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
    bool res=false;
    // if(unif(gen) < exp(-1*delta_E)){
    //     res=true;
    // }
    // cout<<-1*delta_E<<'\n';
    // if(delta_E<0){
    //     res=true;
    // }    
    // if(delta_E!=0){cout<<delta_E<<'\n';}
    
    // coust<<exp(-10*delta_E)<<'\t'<<res<<'\n';
    return (unif(gen) < exp(-0.1*delta_E));//return 1接受 0拒绝
    // return (delta_E<0);//return 1接受 0拒绝
}

//---------------------------boundary---------------------------------

bool check_boundary(Vector3i prop_move1) {
    //find if point in sphere
    if (pow(prop_move1[0],2)+pow(prop_move1[1],2)+pow(prop_move1[2],2)<=pow(radius,2)) { 
        return 1;
    } 
    else {return 0;}
};

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

void tad_move(int thread_num,vector<MatrixXi> &connected_structure_chr,vector<vector<MatrixXi>> &tads_chr,vector<MatrixXf> &tads_center_chr,vector<MatrixXi> &dir_match,int chr_homology,int tad_index,int m){

    //lock
    mutexs[thread_num].lock();

    int tmp;
    int chr_index=chr_homology/2;
    int tad_len=chr_tad_len[chr_index];
    int chr_bin_len=chr_len[chr_index];    
    int tad_mid_len=tad_len/2;

    int start=chr_tad2bound[chr_index][tad_index][0]; 
    int end=chr_tad2bound[chr_index][tad_index][1];

    MatrixXf tads_center_change;  
    MatrixXi change_structure;       
    MatrixXi first_point(1,3);
    MatrixXi end_point(1,3);     
    MatrixXi old_tad=tads_chr[chr_homology][tad_index];  
    MatrixXi tad ;

    int dir=dir_match[chr_homology](0,tad_index);

    uniform_int_distribution<int> uniconnect(0,7);
    int dir_change=uniconnect(gen); 
    while(dir==dir_change){
        dir_change=uniconnect(gen); 
    }

    //update m
    thread_m[thread_num]=m;     

    if(tad_index>=tad_mid_len && tad_index!=tad_len-1){

        //intial dis_matrix_change
        int start_index=homology_tad_start[chr_homology];
        // homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+start_index+tad_index,homology_dis_matrix[thread_num].begin()+start_index+tad_len);        
        homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin(),homology_dis_matrix[thread_num].end());

        first_point=connected_structure_chr[chr_homology](start,all);
        move_connect_tad(first_point,old_tad,tad,dir_change);  
        
        end_point=tad(tad.rows()-1,all);    

        change_structure=connected_structure_chr[chr_homology].bottomRows(chr_bin_len-end);
        end_point-=change_structure(0,all);
        move_tad(end_point,change_structure,1);      

        Matrix<float, 1, 3> tad_center(float(tad.topRows(end-start).col(0).sum())/(end-start),float(tad.topRows(end-start).col(1).sum())/(end-start),float(tad.topRows(end-start).col(2).sum())/(end-start));
        tads_center_change=tads_center_chr[chr_homology].bottomRows(tad_len-tad_index);
        MatrixXf move_point(1,3);
        move_point<<float(end_point(0,0)),float(end_point(0,1)),float(end_point(0,2));  
        move_center(move_point,tads_center_change,1.0);  
        tads_center_change.topRows(1)=tad_center;

        //update homology_dis_matrix_change and caculate energy    
        float energy=get_delta_E(thread_num,chr_homology,tad_index,tads_center_change,tads_center_chr,homology_dis_matrix[thread_num]);
        // cout<<"energy:"<<energy<<'\n';
        // //get sphere_limition
        float sphere_limition=sphere_limit(tads_center_change,alpha)-sphere_limit(tads_center_chr[chr_homology].bottomRows(tad_len-tad_index),alpha);            
        // float sphere_limition=0;
        
        if(accept_move(energy + sphere_limition)){

            //update coordinate               
            connected_structure_chr[chr_homology].bottomRows(chr_bin_len-end)=change_structure; 
            connected_structure_chr[chr_homology].middleRows(start,end-start+1)=tad;                            
        
            //update tad_center   
            tads_center_chr[chr_homology].bottomRows(tad_len-tad_index)=tads_center_change;         

            //update dir_match  
            dir_match[chr_homology](0,tad_index)=dir_change;    

            //update contacts
            update_neighbor_contact(chr_homology,thread_num,m);            

            //update neighbor
            update_dis_matrix(thread_num,chr_homology,tad_index,1);    


        }               
    }  
    else if(tad_index==tad_len-1){

        //intial dis_matrix_change
        int homology_index=homology_tad_start[chr_homology]+tad_index;
        // homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+homology_index,homology_dis_matrix[thread_num].begin()+homology_index+1);
        homology_dis_matrix_change[thread_num]=homology_dis_matrix[thread_num];
        //get tads_center_change
        first_point=connected_structure_chr[chr_homology](start,all);
        move_connect_tad(first_point,old_tad,tad,dir_change);               
        end_point=tad(tad.rows()-1,all);   

        Matrix<float, 1, 3> tad_center(float(tad.topRows(end-start).col(0).sum())/(end-start),float(tad.topRows(end-start).col(1).sum())/(end-start),float(tad.topRows(end-start).col(2).sum())/(end-start));
        tads_center_change=tad_center;    

        float energy=get_delta_E(thread_num,chr_homology,tad_index,tads_center_change,tads_center_chr,homology_dis_matrix[thread_num]);     

        //get sphere_limition
        float sphere_limition=sphere_limit(tads_center_change,alpha)-sphere_limit(tads_center_chr[chr_homology](tad_index,all),alpha);            
        // float sphere_limition=0;
        
        if(accept_move(energy + sphere_limition)){

            //update coordinate 
            connected_structure_chr[chr_homology].middleRows(start,end-start)=tad; 

            //update tad_center 
            tads_center_chr[chr_homology](tad_index,all)=tad_center;            

            //update dir_match  
            dir_match[chr_homology](0,tad_index)=dir_change;    

            //update contacts
            update_neighbor_contact(chr_homology,thread_num,m);               
              
            //update neighbor
            homology_dis_matrix[thread_num]=homology_dis_matrix_change[thread_num];
                    
        }                    

    }      
    else{

        //intial dis_matrix_change
        int start_index=homology_tad_start[chr_homology];
        // homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+start_index,homology_dis_matrix[thread_num].begin()+start_index+tad_index+1);    
        homology_dis_matrix_change[thread_num]=homology_dis_matrix[thread_num];
        move_tad(old_tad.bottomRows(1),old_tad,-1);
        first_point=connected_structure_chr[chr_homology](end,all);

        move_connect_tad(first_point,old_tad,tad,dir_change);
        end_point=tad(0,all);    

        MatrixXi change_structure=connected_structure_chr[chr_homology].topRows(start+1);
        end_point-=change_structure.bottomRows(1);
        move_tad(end_point,change_structure,1);

        tads_center_change=tads_center_chr[chr_homology].topRows(tad_index+1);        
        MatrixXf move_point(1,3);
        move_point<<float(end_point(0,0)),float(end_point(0,1)),float(end_point(0,2));

        move_center(move_point,tads_center_change,1.0); 
        Matrix<float, 1, 3> tad_center(float(tad.topRows(end-start).col(0).sum())/(end-start),float(tad.topRows(end-start).col(1).sum())/(end-start),float(tad.topRows(end-start).col(2).sum())/(end-start));        
        tads_center_change.bottomRows(1)=tad_center;

        float energy=get_delta_E(thread_num,chr_homology,tad_index,tads_center_change,tads_center_chr,homology_dis_matrix[thread_num]);

        //get sphere_limition
        float sphere_limition=sphere_limit(tads_center_change,alpha)-sphere_limit(tads_center_chr[chr_homology].topRows(tad_index+1),alpha);            
        // float sphere_limition=0;
        
        if(accept_move(energy + sphere_limition)){

            //update tad_center   
            tads_center_chr[chr_homology].topRows(tad_index+1)=tads_center_change;             

            // //update coordinate
            connected_structure_chr[chr_homology].middleRows(start,end-start+1)=tad;       
            connected_structure_chr[chr_homology].topRows(start+1)=change_structure;    

            //update dir_match  
            dir_match[chr_homology](0,tad_index)=dir_change;    

            //update contacts
            update_neighbor_contact(chr_homology,thread_num,m);                

            //update neighbor
            update_dis_matrix(thread_num,chr_homology,tad_index,0); 
          
                     
        }
    }

    //unlock
    mutexs[thread_num].unlock();    

    
}

void tad_change(int thread_num,int chr_homology,int tad_index,int m){

    // lock and update change_flags
    tad_change_mutex.lock() ;

    int chr_index=chr_homology/2;
    int tad_len=chr_tad_len[chr_index];
    int tad_mid_len=tad_len/2;
    int chr_bin_len=chr_len[chr_index];

    //initial change thread
    uniform_int_distribution<int> uniIndex(0,tad_len-1);
    tad_index = uniIndex(gen);         

    int change_thread=uniThread(gen);
    while(thread_num==change_thread){
        change_thread=uniThread(gen);
    }   

    //if one of the two thread is running change tad ,break
    if ((tad_change_flags[thread_num]==1) || (tad_change_flags[change_thread]==1)){
        tad_change_mutex.unlock() ;
        return;
    }
    tad_change_flags[thread_num]= 1 ;
    tad_change_flags[change_thread]= 1 ;
    tad_change_mutex.unlock() ;    
  
    (mutexs[change_thread]).lock();

    //update m
    thread_m[thread_num]=m;    

    if(tad_index>=tad_mid_len){

        //intial dis_matrix_change
        int start_index=homology_tad_start[chr_homology];
        // homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+start_index+tad_index,homology_dis_matrix[thread_num].begin()+start_index+tad_len);    
        // homology_dis_matrix_change[change_thread].assign(homology_dis_matrix[change_thread].begin()+start_index+tad_index,homology_dis_matrix[change_thread].begin()+start_index+tad_len);    

        homology_dis_matrix_change[thread_num]=homology_dis_matrix[thread_num];
        homology_dis_matrix_change[change_thread]=homology_dis_matrix[change_thread];

        int start=chr_tad2bound[chr_index][tad_index][0];
        int end=chr_tad2bound[chr_index][tad_len-1][1];    

        //exchange tad
        MatrixXi change1=connected_structure_chr[thread_num][chr_homology].bottomRows(chr_bin_len-start);  
        MatrixXi change2=connected_structure_chr[change_thread][chr_homology].bottomRows(chr_bin_len-start);  
    
        MatrixXi margin1=connected_structure_chr[thread_num][chr_homology](start,all);   
        MatrixXi margin2=connected_structure_chr[change_thread][chr_homology](start,all); 
        MatrixXi margin =margin2-margin1 ;

        move_tad(margin,change1,1);
        move_tad(margin,change2,-1);

        // update tad center         
        MatrixXf tads_center_change1=tads_center_chr[thread_num][chr_homology].bottomRows(tad_len-tad_index);
        MatrixXf tads_center_change2=tads_center_chr[change_thread][chr_homology].bottomRows(tad_len-tad_index);
        
        Matrix<float,1,3> margin_center(float((margin2-margin1)(0,0)),float((margin2-margin1)(0,1)),float((margin2-margin1)(0,2)));

        move_center(margin_center,tads_center_change1,1);
        move_center(margin_center,tads_center_change2,-1);   

        // //get_energy              
        float energy1=get_delta_E(thread_num,chr_homology,tad_index,tads_center_change2,tads_center_chr[thread_num],homology_dis_matrix[thread_num]);
        float energy2=get_delta_E(change_thread,chr_homology,tad_index,tads_center_change1,tads_center_chr[change_thread],homology_dis_matrix[change_thread]);

        //get sphere_limition
        float sphere_limition1=sphere_limit(tads_center_change1,alpha)-sphere_limit(tads_center_chr[change_thread][chr_homology].bottomRows(tad_len-tad_index),alpha);            
        float sphere_limition2=sphere_limit(tads_center_change2,alpha)-sphere_limit(tads_center_chr[thread_num][chr_homology].bottomRows(tad_len-tad_index),alpha);    

        if(accept_move(energy1+energy2+sphere_limition1+sphere_limition2)){
            // cout<<thread_num<<'\t'<<change_thread<<'\t'<<tad_index<<'\n';

            //update tad_center
            tads_center_chr[thread_num][chr_homology].bottomRows(tad_len-tad_index)=tads_center_change2;
            tads_center_chr[change_thread][chr_homology].bottomRows(tad_len-tad_index)=tads_center_change1;

            //update location
            connected_structure_chr[thread_num][chr_homology].bottomRows(chr_bin_len-start)=change2;  
            connected_structure_chr[change_thread][chr_homology].bottomRows(chr_bin_len-start)=change1;          ;

            //update contacts
            update_neighbor_contact(chr_homology,thread_num,m);     
            update_neighbor_contact(chr_homology,change_thread,m);    

            //update neighbor
            update_dis_matrix(thread_num,chr_homology,tad_index,1); 
            update_dis_matrix(change_thread,chr_homology,tad_index,1); 
             
        }
    }

    else{

        //intial dis_matrix_change
        int start_index=homology_tad_start[chr_homology];
        // homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+start_index,homology_dis_matrix[thread_num].begin()+start_index+tad_index+1);    
        // homology_dis_matrix_change[change_thread].assign(homology_dis_matrix[change_thread].begin()+start_index,homology_dis_matrix[change_thread].begin()+start_index+tad_index+1);  

        homology_dis_matrix_change[thread_num]=homology_dis_matrix[thread_num];
        homology_dis_matrix_change[change_thread]=homology_dis_matrix[change_thread];

        int start=0;
        int end=chr_tad2bound[chr_index][tad_index][1];   

        //exchange tad
        MatrixXi change1=connected_structure_chr[thread_num][chr_homology].topRows(end+1);  
        MatrixXi change2=connected_structure_chr[change_thread][chr_homology].topRows(end+1);  
    
        MatrixXi margin1=connected_structure_chr[thread_num][chr_homology](end,all);   
        MatrixXi margin2=connected_structure_chr[change_thread][chr_homology](end,all); 
        MatrixXi margin =margin2-margin1;

        move_tad(margin,change1,1);
        move_tad(margin,change2,-1);

        // update tad center         
        MatrixXf tads_center_change1=tads_center_chr[thread_num][chr_homology].topRows(tad_index+1);
        MatrixXf tads_center_change2=tads_center_chr[change_thread][chr_homology].topRows(tad_index+1);
        Matrix<float,1,3> margin_center(float((margin)(0,0)),float((margin)(0,1)),float((margin)(0,2)));

        move_center(margin_center,tads_center_change1,1);
        move_center(margin_center,tads_center_change2,-1);

        float energy1=get_delta_E(thread_num,chr_homology,tad_index,tads_center_change2,tads_center_chr[thread_num],homology_dis_matrix[thread_num]);
        float energy2=get_delta_E(change_thread,chr_homology,tad_index,tads_center_change1,tads_center_chr[change_thread],homology_dis_matrix[change_thread]);

        float sphere_limition1=sphere_limit(tads_center_change1,alpha)-sphere_limit(tads_center_chr[change_thread][chr_homology].topRows(tad_index+1),alpha);            
        float sphere_limition2=sphere_limit(tads_center_change2,alpha)-sphere_limit(tads_center_chr[thread_num][chr_homology].topRows(tad_index+1),alpha);    

        if(accept_move(energy1+energy2+sphere_limition1+sphere_limition2)){        

            // update tad_center
            tads_center_chr[thread_num][chr_homology].topRows(tad_index+1)=tads_center_change2;
            tads_center_chr[change_thread][chr_homology].topRows(tad_index+1)=tads_center_change1;

            //update location
            connected_structure_chr[thread_num][chr_homology].topRows(end+1)=change2;  
            connected_structure_chr[change_thread][chr_homology].topRows(end+1)=change1;  

            //update contacts
            update_neighbor_contact(chr_homology,thread_num,m);     
            update_neighbor_contact(chr_homology,change_thread,m);                 
        
            //update neighbor
            update_dis_matrix(thread_num,chr_homology,tad_index,0); 
            update_dis_matrix(change_thread,chr_homology,tad_index,0); 
         
        
        }       
    }    

    //unlock
    (mutexs[change_thread]).unlock();   

    // update tad_match
    for(int i=tad_index;i<tad_len;++i){
        // cout<<tad_match_chr[chr_homology](thread_num,i)<<'\t'<<tad_match_chr[chr_homology](change_thread,i)<<'\n';
        int tmp=tad_match_chr[chr_homology](thread_num,i);
        tad_match_chr[chr_homology](thread_num,i)=tad_match_chr[chr_homology](change_thread,i);
        tad_match_chr[chr_homology](change_thread,i)=tmp;
    }

    // update change_flags
    tad_change_mutex.lock();
    tad_change_flags[thread_num]= 0 ;
    tad_change_flags[change_thread]= 0 ;
    tad_change_mutex.unlock();     
}

#endif //INVERSE_3D_NEW_MOVES_H