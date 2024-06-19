#ifndef INVERSE_3D_MC_H//Inverse_Monte_Carlo
#define INVERSE_3D_MC_H


#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "Moves_liner.h"
#include "global.h"
#include "initial_polymer.h"


void move(std::vector<vector<int>> &polymer,std::vector<int> &site_index2position, std::vector<int> &site_available,int thread_num, int m,int update,int site){ //performs a single Monte Carlo step
    //add and remove have low frequence
    
    float frequence=unif(gen);

    if (frequence<=0.4){
        kink_move(polymer,site_index2position, site_available , site, thread_num,m,update);
    }
    else if (frequence<=0.8){
        loop_move(polymer,site_index2position, site_available , site, thread_num,m,update);
    }
    else if (frequence<=0.9){
       loop_reduce_move(polymer,site_index2position, site_available , site, thread_num,m,update);
    }
    else{
       loop_add_move(polymer,site_index2position, site_available , site, thread_num,m,update);
    }    
}

void run_burnin(int thread_num, int mc_moves) { //burns in the polymer configurations
    std::uniform_int_distribution<int> unisite(0,pol_length-1);
    
    for (int m = 1; m < mc_moves; m++) {
        int site = unisite(gen);
        move(polymer[thread_num],site_index2position[thread_num], site_available[thread_num], thread_num, m, 0,site);
    }
}

void run(int thread_num, int mc_moves) {

    //reset contact frequencies before starting new forward round
    unordered_map<pair<int, int>, float, pair_hash>().swap(one_tad_contact[thread_num]);
    unordered_map<int, vector<int>>().swap(tad_neighbor_contact[thread_num]);
    unordered_map<pair<int, int>, float, pair_hash>().swap(tad_total_contacts_neighbor[thread_num]); 

    if(thread_num==0){
        cout<<neighbor_contact_data[thread_num].size()<<'\n';
    }
    

    //将m更新为0，但是接触不要变  
    for(auto &hash:neighbor_contact_data[thread_num]){
        hash.second.m=0;
    }
    for(auto &hash:one_tad_contact[thread_num]){
        hash.second=0;
    }    

       
    tad_total_contacts_s[thread_num]=tad_total_contacts_s[thread_num]*0;
    // neighbor_contact_m_mat[l]=Eigen::MatrixXi::Zero(tad_bin_num,tad_bin_num);           

    std::uniform_int_distribution<int> unisite(0,pol_length-1);
    
    for (int m = 1; m < mc_moves; m++) {    //performs a forward polymer simulation
        int site = unisite(gen);
        move(polymer[thread_num],site_index2position[thread_num], site_available[thread_num], thread_num, m,  1,site);
    }

    //add_neighbor_contact
    typedef Eigen::Triplet<float> T;
    std::vector<T> c;
    for(auto &hash:tad_total_contacts_neighbor[thread_num]){
        c.push_back(T(hash.first.first,hash.first.second,hash.second));
    }

    Eigen::SparseMatrix<float> mat(tad_bin_num,tad_bin_num);
    mat.setFromTriplets(c.begin(),c.end()); 
    tad_total_contacts_s[thread_num]+=mat; 

    //add remain neighbor contact
    std::vector<T> c2;
    float remain_contact;
    float dis;
    int sitei,sitej,m_pre; 
    for(auto &hash:neighbor_contact_data[thread_num]){
        sitei=hash.first.first/6;
        sitej=hash.first.second/6;
        if(abs(sitei-sitej)<2){continue;}
        m_pre=hash.second.m;
        dis=hash.second.dis;
        if(dis==0){continue;}
        remain_contact=(mc_moves-m_pre)/dis/dis/dis;

        if(mc_moves-m_pre<0){cout<<"error11:"<<m_pre<<'\n';}
        
        // if(abs(sitei-sitej)>=200){cout<<sitei<<'\t'<<sitej<<'\t'<<m_pre<<'\t'<<dis<<'\t'<<remain_contact<<'\n';}
        c2.push_back(T(sitei,sitej,remain_contact));
    }

    Eigen::SparseMatrix<float> mat2(tad_bin_num,tad_bin_num);
    mat2.setFromTriplets(c2.begin(),c2.end());  
    tad_total_contacts_s[thread_num]+=mat2;     

    //add remain bin contact
    std::vector<T> c3;
    for(auto &hash:one_tad_contact[thread_num]){
        sitei=hash.first.first/6;
        sitej=hash.first.second/6;
        if(abs(sitei-sitej)<2){continue;}
        m_pre=hash.second;
        remain_contact=2*(mc_moves-m_pre);
        c3.push_back(T(sitei,sitej,remain_contact));
    }   

    Eigen::SparseMatrix<float> mat3(tad_bin_num,tad_bin_num);
    mat3.setFromTriplets(c3.begin(),c3.end());  
    tad_total_contacts_s[thread_num]+=mat3;        

}

#endif