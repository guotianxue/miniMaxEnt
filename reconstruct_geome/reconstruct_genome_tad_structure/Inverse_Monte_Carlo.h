#ifndef INVERSE_3D_MC_H//Inverse_Monte_Carlo
#define INVERSE_3D_MC_H


#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include "Moves_liner.h"
#include "Chrom_Moves.h"
#include "global.h"
#include "initial_polymer.h"

//===============reconstruct structures============
void tads_move(int thread_num,int m){ //performs a single Monte Carlo step

    float frequence=unif(gen);

    int homology_tad_num=homology_tad2chr_tad.size();
    uniform_int_distribution<int> unitad(0,homology_tad_num-1);
    int homology_tad_index=unitad(gen);
    int chr_homology=homology_tad2chr_tad[homology_tad_index][0];
    int tad_index=homology_tad2chr_tad[homology_tad_index][1];

    // tad_move(thread_num,connected_structure_chr[thread_num],tads_chr[thread_num],tads_center_chr[thread_num],dir_match[thread_num],chr_homology,tad_index,m);
    
    if (frequence<0.9){
        tad_move(thread_num,connected_structure_chr[thread_num],tads_chr[thread_num],tads_center_chr[thread_num],dir_match[thread_num],chr_homology,tad_index,m);
    }
    else{   
        // tad_change(thread_num,chr_homology,tad_index,m);
    }

}

void chrom_tads_move(int thread_num,int m){ //performs a single Monte Carlo step

    uniform_int_distribution<int> unimove(0,1);
    int move_type=unimove(gen);

    uniform_int_distribution<int> unichrom(0,chr_homology_num-1);
    int chr_homology=unichrom(gen);
    
    // int homology_tad_num=homology_tad2chr_tad.size();
    // uniform_int_distribution<int> unitad(0,homology_tad_num-1);
    // int homology_tad_index=unitad(gen);
    // int chr_homology=homology_tad2chr_tad[homology_tad_index][0];
    // int tad_index=homology_tad2chr_tad[homology_tad_index][1];   

    // chromosome_move(thread_num,connected_structure_chr[thread_num],tads_chr[thread_num],tads_center_chr[thread_num],dir_match[thread_num],chr_homology,m); 

    // chromosome_rotate(thread_num,connected_structure_chr[thread_num],tads_chr[thread_num],tads_center_chr[thread_num],dir_match[thread_num],chr_homology,m);    

    if(move_type==0){
        chromosome_move(thread_num,connected_structure_chr[thread_num],tads_chr[thread_num],tads_center_chr[thread_num],dir_match[thread_num],chr_homology,m);
    }
    else{
        chromosome_rotate(thread_num,connected_structure_chr[thread_num],tads_chr[thread_num],tads_center_chr[thread_num],dir_match[thread_num],chr_homology,m);
    }

}

void run_tads(int thread_num, int mc_moves) { //burns in the polymer configurations

    for (int m = 1; m < mc_moves; m++) {
        tads_move(thread_num,m);
    
        // chrom_tads_move(thread_num,m);
    } 

    // read in remaining contacts at the end of forward simulation

    for(int k=0;k<chr_homology_num;k++){    
        update_neighbor_contact(k,thread_num,mc_moves-1);          
    }    
 

    // //dis_matrix
    // for(int k=0;k<tads_num;k++){
    //     for(int j=0;j<tads_num;j++){

    //         if(homology_contact_matrix[0][k][j]<0){
    //             cout<<k<<'\t'<<j<<'\t'<<homology_dis_matrix[0][k][j]<<'\n';
    //         }
    //     }
    // }   

    // //输出homology_contact_matrix看一下结果
    // string dis_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure/dis.txt";
    // ofstream distance(dis_file);
    // for(int i=0;i<500;i++){//2958
    //     for(int j=0;j<500;j++){
    //         distance<<homology_dis_matrix[0][i][j]<<'\t';
    //         // distance<<homology_dis_matrix[0][i][j]<<'\t';
    //     }
    //     distance<<'\n';
    // }
    // distance.close();               


    // if(thread_num==0){  

    //     cout<<"##########################"<<'\n';
    //     cout<<tads_center_chr[thread_num][43]<<'\n';
    //     cout<<'\n';
    //     cout<<'\n';    
      
    // }    

}


void run_chrom_tads(int thread_num, int mc_moves) { //burns in the polymer configurations

    for (int m = 1; m < mc_moves; m++) {
        // tads_move(thread_num,m);
    
        chrom_tads_move(thread_num,m);
    } 

    // read in remaining contacts at the end of forward simulation

    for(int k=0;k<chr_homology_num;k++){    
        update_chrom_neighbor_contact(k,thread_num,mc_moves);          
    }    
 
}

#endif