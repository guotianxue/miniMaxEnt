//
// Created by Joris on 09/07/2018.
//
#ifndef INVERSE_3D_NEW_MOVES_H
#define INVERSE_3D_NEW_MOVES_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <queue>
#include <stack>
#include "Energy_changes_liner.h"

using namespace std;

//test
int reduce_count=0,add_count=0,kink_count=0,loop_count=0,skip_count=0;

uniform_int_distribution<int> unidir(0,2);
uniform_int_distribution<int> unidir_loop(1,5);

//---------------------------accept move---------------------------------

bool accept_move(float delta_E){//如果能量变小了，接受Move。如果能量变大了，以一定概率接受Move
    float tmp=unif(gen);
    float tmp2=exp(-delta_E);
    delta_E=delta_E*100;
    return ((delta_E <=0)||(unif(gen) < exp(-delta_E)));//return 1接受 0拒绝
}

//---------------------------boundary---------------------------------

bool check_boundary(vector<int> prop_move1) {
    //find if point in sphere
    if (pow(prop_move1[0],2)+pow(prop_move1[1],2)+pow(prop_move1[2],2)<=pow(radius,2)) { 
        return 1;
    } 
    else {return 0;}
};

float sphere_limit(vector<int> move,float alpha){
    float sphere_out=0;

    float distance=sqrt(pow(move[0],2)+pow(move[1],2)+pow(move[2],2));

    if (distance>radius){
        sphere_out+=(distance-radius)*alpha;
    }

    return sphere_out;
}

//
bool is_equal_loc(vector<int> loc1,vector<int> loc2){
    int len1=loc1.size();
    int len2=loc2.size();
    if(len1!=len2){return false;}
    for(int i=0;i<len1;i++){
        if(loc1[i]!=loc2[i]){
            return false;
        }
    }
    return true;
}

vector<int> add_loc(vector<int> loc1,vector<int> loc2,int k1,int k2){
    int len1=loc1.size();
    int len2=loc2.size();
    vector<int> res(len1);

    for(int i=0;i<len1;i++){
        res[i]=k1*loc1[i]+k2*loc2[i];
    }
    return res;
}

//-------------------apply to the condition which the two site is far----------------------
void deal_neighbor( int number_of_threads,int site,int pol_length,int m,vector<int> prop_move){
    
    //initial tad_neighbor_contact_new
    if (site%6!=0){return;}

    int x=prop_move[0];int y=prop_move[1];int z=prop_move[2];
    find_neighbor(x,y,z,site,m,polymer[number_of_threads],polymer_bin[number_of_threads],site_index2position[number_of_threads],tad_neighbor_contact_new[number_of_threads],neighbor_contact_dis[number_of_threads]);

    return ;

}


//------------------------------------four moves---------------------------------------------------
void boundary_move(vector<vector<int>> &polymer, vector<int> &site_index2position, vector<int> &site_available, int site,int thread_num, int m){//还没加update,过一会加
    vector<int> prop_move1(3);
    int flag;
    if (site==pol_length-1){
        flag=1;
    }
    else if  (site==0){
        flag=-1;
    }    
    int direction = unidir_loop(gen);//1-5
    vector<int> rotate_vector(3);
    for (int i = 0; i<3;i++){ // rotate loop in one of 5 possible new directions
        int k=(i+direction/2)%3;
        rotate_vector[k] = (-2*(direction%2)+1)*(polymer[site_index2position[site]][i]-polymer[site_index2position[site]-flag][i]);
        prop_move1[k]=polymer[site_index2position[site]-flag][k]+rotate_vector[k]; 
    }        

    // deal with neighbor contact           
    deal_neighbor(thread_num,site,pol_length,m,prop_move1);     

    float energy_change=delta_E_other(polymer, site_index2position,site_available, site, pol_length, prop_move1,thread_num);
    // float distace_margin=sphere_limit(prop_move1,1)-sphere_limit(polymer[site_index2position[site]],1);

    // if((accept_move(energy_change)==1 && distace_margin<=0) || distace_margin<0){
    //  if(accept_move(energy_change)==1 && check_nucleus(prop_move1[0],prop_move1[1],prop_move1[2]) ){ 
    if(accept_move(energy_change)==1 && check_nucleus(prop_move1[0],prop_move1[1],prop_move1[2]) ){ 
        //update contact
        update_contact_location(polymer, site, pol_length, prop_move1,polymer[site_index2position[site]],thread_num,m);
           
        update_contact_location_neighbor(site,thread_num,m);    
        
        //update polymer
        polymer[site_index2position[site]]=prop_move1; 
        if(site%6==0){
            polymer_bin[thread_num][site/6]=prop_move1;
        }  
    }

}

void kink_move(vector<vector<int>> &polymer, vector<int> &site_index2position, vector<int> &site_available, int site, int thread_num, int m,int update){
    
    vector<int> site_use;
    int flag=0;
    while(site<pol_length-1 and flag<3){
        if (site_available[site]==1){
            site_use.emplace_back(site);
            ++flag;            
        }
        site+=1;
    }

    if (site==pol_length-1 || site==0){
        boundary_move(polymer,site_index2position, site_available,site,thread_num,m);
    }
    else if (site_use.size()==3 and site_use[2]<pol_length &&  (is_equal_loc(polymer[site_index2position[site_use[2]]],polymer[site_index2position[site_use[0]]])==false)  &&  (is_equal_loc(polymer[site_index2position[site_use[2]]] , add_loc(polymer[site_index2position[site_use[1]]],polymer[site_index2position[site_use[0]]],2,-1))==false)){
        //第一个条件：两个位点后的坐标不和原位点重合    第二个约束条件：确保这三个位点不在一条直线上，得有转角
        kink_count+=1;
        vector<int> prop_move1(3);
        for(int i=0;i<3;i++){
            prop_move1[i] = polymer[site_index2position[site_use[0]]][i] + polymer[site_index2position[site_use[2]]][i] - polymer[site_index2position[site_use[1]]][i];//site+1位置的kink_move之后的坐标
        }
        

        //deal with neighbor contact           
        deal_neighbor(thread_num,site_use[1],pol_length,m,prop_move1);              

        // float distace_margin=sphere_limit(prop_move1,1)-sphere_limit(polymer[site_index2position[site_use[1]]],1);      
        float energy_change=delta_E_other(polymer, site_index2position ,site_available,site_use[1], pol_length, prop_move1,thread_num);

        // if((accept_move(energy_change)==1 && distace_margin<=0) || distace_margin<0){        
        //  if((accept_move(energy_change)==1 && check_nucleus(prop_move1[0],prop_move1[1],prop_move1[2]) )){  
        if((accept_move(energy_change)==1 && check_nucleus(prop_move1[0],prop_move1[1],prop_move1[2]) )){  
            //update contact
            update_contact_location(polymer, site_use[1], pol_length, prop_move1,polymer[site_index2position[site_use[1]]],thread_num,m);
               
            update_contact_location_neighbor(site_use[1],thread_num,m);    
            
            //update polymer
            polymer[site_index2position[site_use[1]]]=prop_move1;      
            if(site%6==0){
                polymer_bin[thread_num][site/6]=prop_move1;
            }                      
        }      

    }       
}


void loop_move(vector<vector<int>> &polymer, vector<int> &site_index2position, vector<int> &site_available, int site, int thread_num,int m,int update){
    
    vector<int> site_use;
    int flag=0;
    while(site<pol_length-1 and flag<3){
        if (site_available[site]==1){
            site_use.emplace_back(site);
            ++flag;
        }
        site+=1;
    }    
    if (site==pol_length-1 || site==0){
        boundary_move(polymer,site_index2position, site_available,site,thread_num,m);

    }
    else if (site_use.size()==3 && site_use[2]<pol_length && is_equal_loc(polymer[site_index2position[site_use[2]]],polymer[site_index2position[site_use[0]]])==true){
        loop_count+=1;
        int direction = unidir_loop(gen);
        vector<int> rotated_vector(3); vector<int> prop_move1(3);

        for (int i = 0; i<3;i++){ // rotate loop in one of 5 possible new directions
            int k=(i+direction/2)%3;
            rotated_vector[k] = (-2*(direction%2)+1)*(polymer[site_index2position[site_use[1]]][i]-polymer[site_index2position[site_use[0]]][i]);
            prop_move1[k] = polymer[site_index2position[site_use[0]]][k] + rotated_vector[k];
        }
        

        //deal with neighbor contact           
        deal_neighbor(thread_num,site_use[1],pol_length,m,prop_move1);    

        float distace_margin=sphere_limit(prop_move1,1)-sphere_limit(polymer[site_index2position[site_use[1]]],1);      
        float energy_change=delta_E_other(polymer, site_index2position ,site_available,site_use[1], pol_length, prop_move1,thread_num);
        // if((accept_move(energy_change)==1 && distace_margin<=0) || distace_margin<0){
        //  if((accept_move(energy_change)==1 && check_nucleus(prop_move1[0],prop_move1[1],prop_move1[2]) )){ 
        if((accept_move(energy_change)==1) && check_nucleus(prop_move1[0],prop_move1[1],prop_move1[2])){
            //update contact
            update_contact_location(polymer, site_use[1], pol_length, prop_move1,polymer[site_index2position[site_use[1]]],thread_num,m);
            
            update_contact_location_neighbor(site_use[1],thread_num,m);    
            
            //update polymer
            polymer[site_index2position[site_use[1]]]=prop_move1;  
            if(site%6==0){
                polymer_bin[thread_num][site/6]=prop_move1;
            }              
        }
        else{
            skip_count+=1;
        } 
    }
}


void loop_reduce_move(vector<vector<int>> &polymer,vector<int> &site_index2position, vector<int> &site_available, int site, int thread_num,int m,int update){

    vector<int> site_use;
    int flag=1;
    int direction = leftOright(gen);    
    if (direction){flag=-1;}
    if (site!=0 && site%6==0){
        site_use.emplace_back(site);
        if (site_available[site+flag]) {
            site_use.emplace_back(site+flag);
        }
        if (site_available[site+flag*2]) {
            site_use.emplace_back(site+flag*2);
        }     
    }
         
    if (site==pol_length-1 || site==0 ){
        boundary_move(polymer,site_index2position, site_available,site,thread_num,m);
    }
    else if ( site_use.size()==3 && is_equal_loc(polymer[site_index2position[site_use[2]]] , polymer[site_index2position[site_use[0]]])==false){     
        vector<int> margin(3);
        for(int i=0;i<3;i++){
            margin[i] = polymer[site_index2position[site_use[2]]][i]- polymer[site_index2position[site_use[0]]-flag][i] ;
        }
         
        if(abs(margin[0])+abs(margin[1])+abs(abs(margin[2]))==1){
            vector<int> prop_move=polymer[site_index2position[site_use[2]]];
            
            //deal with neighbor contact
            deal_neighbor(thread_num,site_use[0],pol_length,m,prop_move);              

            float distace_margin=sphere_limit(prop_move,1)-sphere_limit(polymer[site_index2position[site_use[0]]],1);     

            float energy_change=delta_E_other(polymer, site_index2position,site_available, site_use[0], pol_length,prop_move,thread_num);
            // if((energy_change<0 && distace_margin<=0) || distace_margin<0){//if(accept_move(energy_change)==1){
            // if((energy_change<0 && check_nucleus(prop_move[0],prop_move[1],prop_move[2]) )){   
            if(energy_change<0 && check_nucleus(prop_move[0],prop_move[1],prop_move[2])){
                ++reduce_count;
                //update contact and location
                update_contact_location(polymer, site_use[0], pol_length,prop_move,polymer[site_index2position[site_use[0]]],thread_num,m);
                polymer[site_index2position[site_use[0]]]=prop_move;
                polymer_bin[thread_num][site_use[0]/6]=prop_move;
  
                update_contact_location_neighbor(site_use[0],thread_num,m); 

                //update site_index2position and site_available and position2site_index
                int pre_site1=site_index2position[site_use[1]];
                int pre_site2=site_index2position[site_use[2]];

                site_available[site_use[1]]=0;
                site_available[site_use[2]]=0;      
                    
                if (flag==1){
                    site_index2position[site_use[1]]-=1; 
                    for (int i=site_use[2];i<pol_length;++i){
                        site_index2position[i]-=2;
                    }                  
                }
                else{
                    site_index2position[site_use[2]]-=1; 
                    for (int i=site_use[1];i<pol_length;++i){
                        site_index2position[i]-=2;
                    }                          
                }

                //update polymer 
                if (flag>0){polymer.erase(polymer.begin()+pre_site1, polymer.begin()+pre_site2+1);}
                else{polymer.erase(polymer.begin()+pre_site2, polymer.begin()+pre_site1+1);}         
            
            }
        }
    }    
}


void loop_add_move(vector<vector<int>> &polymer, vector<int> &site_index2position, vector<int> &site_available, int site, int thread_num,int m,int update){

    vector<int> site_use;
    int flag=1;
    int direction = leftOright(gen);    
    if (direction){flag=-1;}
    if (site!=0 && site%6==0){
        site_use.emplace_back(site);
        if (site_available[site+flag]==0) {
            site_use.emplace_back(site+flag);
        }
        if (site_available[site+flag*2]==0) {
            site_use.emplace_back(site+flag*2);
        }     
    }  

    if (site==pol_length-1 || site==0 ){//or site==0
        boundary_move(polymer,site_index2position, site_available,site,thread_num,m);
    }
    else if (site_use.size()==3){
        
        vector<int> pre_site;
        vector<int> now_site;
        vector<int> prop_move(3);
        if (site_available[site_use[0]-flag]==0){ pre_site=polymer[site_index2position[site_use[0]-3*flag]]; }
        else{pre_site=polymer[site_index2position[site_use[0]-flag]];}
        now_site=polymer[site_index2position[site_use[0]]];
        //确定loop的方向
        int direction_loop = unidir_loop(gen);
        int move_index=0;
        int move_num=0;
        vector<int> rotate_vector(3);


        for (int i = 0; i<3;i++){ // rotate loop in one of 5 possible new directions
            int k=(i+direction_loop/2)%3;
            rotate_vector[k] = (-2*(direction_loop%2)+1)*(now_site[i]-pre_site[i]);
            prop_move[k]=pre_site[k]+rotate_vector[k]; 
        }        

        //deal with neighbor contact
        deal_neighbor(thread_num,site_use[0],pol_length,m,prop_move);  

        float distace_margin=sphere_limit(prop_move,1)-sphere_limit(polymer[site_index2position[site_use[0]]],1);     
        float energy_change=delta_E_other(polymer, site_index2position,site_available,site_use[0], pol_length,prop_move,thread_num);

        // if((energy_change<0 && distace_margin<=0) || distace_margin<0){         
        // if((energy_change<0 && check_nucleus(prop_move[0],prop_move[1],prop_move[2]) )){   
        if(energy_change<0 && check_nucleus(prop_move[0],prop_move[1],prop_move[2]) ){  
            add_count++;
            //update contact and location
            update_contact_location(polymer, site_use[0], pol_length,prop_move,polymer[site_index2position[site_use[0]]],thread_num,m);
  
            update_contact_location_neighbor(site_use[0],thread_num,m); 

            //update polymer
            polymer[site_index2position[site_use[0]]]=prop_move;      
            polymer_bin[thread_num][site_use[0]/6]=prop_move;      

            //update site_index2position and site_available and position2site_index
            site_available[site_use[1]]=1;
            site_available[site_use[2]]=1;    

            //flag direction
            if (flag==1){
                site_index2position[site_use[1]]+=1; 
                for (int i=site_use[2];i<pol_length;++i){
                    site_index2position[i]+=2;
                }                  
            }
            else{
                site_index2position[site_use[2]]+=1; 
                for (int i=site_use[1];i<pol_length;++i){
                    site_index2position[i]+=2;
                }                          
            }            
            //update polymer 
            //direction
            if (flag==1){
                polymer.insert(polymer.begin()+site_index2position[site_use[1]], pre_site); 
                polymer.insert(polymer.begin()+site_index2position[site_use[2]], now_site);   
            }
            else{
                polymer.insert(polymer.begin()+site_index2position[site_use[2]], now_site); 
                polymer.insert(polymer.begin()+site_index2position[site_use[1]], pre_site);  
            }
        }
                      
    }
          
}


#endif //INVERSE_3D_NEW_MOVES_H