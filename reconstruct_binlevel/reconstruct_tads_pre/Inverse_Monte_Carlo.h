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
    int action = unimove(gen);

    if (action==0){
        kink_move(polymer,site_index2position, site_available , site, thread_num,m,update);
    }
    else if (action==1){
        loop_move(polymer,site_index2position, site_available , site, thread_num,m,update);
    }
    else if (action==2){
       loop_reduce_move(polymer,site_index2position, site_available , site, thread_num,m,update);
    }
    else if (action==3){
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

}

#endif