#ifndef INVERSE_3D_MC_H//Inverse_Monte_Carlo
#define INVERSE_3D_MC_H


#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "Moves_liner.h"
#include "global.h"
#include "initial_polymer.h"

//===============reconstruct structures============
void tads_move(int thread_num,int m){ //performs a single Monte Carlo step
    // float frequence=unif(gen);
    float action = unimove(gen)/10.0;

    if (action<0.5){
        tad_translocation(thread_num,connected_structure[thread_num],is_move_active[thread_num],m);
    }
    else if (action>=0.5 && action<0.8){
        tad_rotate(thread_num,connected_structure[thread_num],is_move_active[thread_num],m);
    }
    else {
        tad_change(thread_num,m);
    }  
}

void run_tads(int thread_num, int mc_moves) { //burns in the polymer configurations

    for (int m = 0; m < mc_moves; m++) {
        tads_move(thread_num,m);
    } 

}

#endif