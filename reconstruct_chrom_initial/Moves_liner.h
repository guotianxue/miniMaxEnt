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

bool accept_move(double delta_E){//如果能量变小了，接受Move。如果能量变大了，以一定概率接受Move
    double tmp=unif(gen);
    double tmp2=exp(-delta_E);
    return (delta_E <=0);
    // return ((delta_E <=0)||(unif(gen) < exp(-delta_E)));//return 1接受 0拒绝
}

//---------------------------boundary---------------------------------

bool check_boundary(vector<int> prop_move) {
    //find if point in sphere
    if (pow(prop_move[0],2)+pow(prop_move[1],2)+pow(prop_move[2],2)<=pow(radius,2)) { 
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



//------------------------------------four moves---------------------------------------------------
void chrom_move(vector<vector<int>> &polymer,  int thread_num){
    
    std::uniform_int_distribution<int> unidir(0,5);
    std::uniform_int_distribution<int> unichro(0,chr_num-1);
    std::uniform_int_distribution<int> unihomo(0,1);
    vector<vector<int>> direction={{0,0,1},{0,0,-1},{0,1,0},{0,-1,0},{1,0,0},{-1,0,0}};

    int dir_index=unidir(gen);
    int chromosome_index=unichro(gen);
    int homology=unihomo(gen);
    int index=chromosome_index*2+homology;
    vector<int> pre_move=polymer[index];
    vector<int> prop_move=add_loc(polymer[index],direction[dir_index],1,1);     
 
    float energy_change=delta_E_chrom(polymer, pre_move, prop_move,chromosome_index,homology);

    if(check_boundary(prop_move) && accept_move(energy_change)==1){               
        //update polymer
        polymer[chromosome_index*2+homology]=prop_move;                        
    }    

    update_chrom_contact(polymer,thread_num);  
}



#endif //INVERSE_3D_NEW_MOVES_H