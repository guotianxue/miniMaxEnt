#ifndef INVERSE_3D_MC_H//Inverse_Monte_Carlo
#define INVERSE_3D_MC_H


#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "Moves_liner.h"
#include "global.h"
#include "initial_polymer.h"

void run(int m,int thread_num) {
    for (int i=0;i<m;i++){
        chrom_move(polymer_chrom[thread_num], thread_num);
    }
   
}

#endif