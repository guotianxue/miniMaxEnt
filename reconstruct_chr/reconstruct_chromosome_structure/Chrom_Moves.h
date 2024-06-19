#ifndef Chrom_MOVES_H
#define Chrom_MOVES_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>
#include <mutex>
#include <chrono>
#include "Energy_changes_liner.h"
#include "Moves_liner.h"
#include "initial_polymer.h"
#include "Chrom_Energy_change.h"


using namespace std;
using namespace Eigen;

//============chromosome level============

void chrom_move(Eigen::MatrixXf &chrom,int move_index){
    
    vector<Eigen::Matrix<float,1,3>> translocation={{0,0,1},{0,0,-1},{0,1,0},{0,-1,0},{1,0,0},{-1,0,0}}; 
    Eigen::MatrixXf move_point=translocation[move_index];
    
    Eigen::MatrixXf margin=move_point.replicate(chrom.rows(), 1);
    chrom=chrom + margin;                
}

void chrom_move(Eigen::MatrixXi &chrom,int move_index){
    
    vector<Eigen::Matrix<int,1,3>> translocation={{0,0,1},{0,0,-1},{0,1,0},{0,-1,0},{1,0,0},{-1,0,0}}; 
    Eigen::MatrixXi move_point=translocation[move_index];

    Eigen::MatrixXi margin=move_point.replicate(chrom.rows(), 1);
    chrom=chrom + margin;                    
}

void chrom_rotate(Eigen::MatrixXf &chrom,int rotate_index){
    vector<Eigen::Matrix<float,1,3>> rotate={{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};    
    Eigen::MatrixXf rotate_mat;
    //rotate
    rotate_mat=rotate[rotate_index].replicate(chrom.rows(), 1);
    chrom=chrom.cwiseProduct(rotate_mat);   //dot                  
}

void chrom_rotate(Eigen::MatrixXi &chrom,int rotate_index){
    vector<Eigen::Matrix<int,1,3>> rotate={{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};    
    Eigen::MatrixXi rotate_mat;
    //rotate
    rotate_mat=rotate[rotate_index].replicate(chrom.rows(), 1);
    chrom=chrom.cwiseProduct(rotate_mat);   //dot                  
}

void rotate_center(Eigen::MatrixXf move_point,Eigen::MatrixXf &tads,int dir){
    vector<Eigen::Matrix<float,1,3>> rotate={{1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};    
    Eigen::MatrixXf margin=move_point.replicate(tads.rows(), 1);
    Eigen::MatrixXf rotate_mat;
    //rotate
    rotate_mat=rotate[dir].replicate(tads.rows(), 1);
    tads=(tads-margin).cwiseProduct(rotate_mat) +margin;   //dot   
}

void rotate_center(Eigen::MatrixXi move_point,Eigen::MatrixXi &tads,int dir){
    vector<Eigen::Matrix<int,1,3>> rotate={{1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};    
    Eigen::MatrixXi margin=move_point.replicate(tads.rows(), 1);
    Eigen::MatrixXi rotate_mat;
    //rotate
    rotate_mat=rotate[dir].replicate(tads.rows(), 1);
    tads=(tads-margin).cwiseProduct(rotate_mat) +margin;   //dot   
}

//=======================chromosome level move======================

void chromosome_move(int thread_num,vector<MatrixXi> &connected_structure_chr,vector<vector<MatrixXi>> &tads_chr,vector<MatrixXf> &tads_center_chr,vector<MatrixXi> &dir_match,int chr_index,int m){

    uniform_int_distribution<int> unimove(0,5); 
    int move_index=unimove(gen);
    int tad_len=chr_tad_len[chr_index];

    int chr_homology1=chr_index*2;
    int chr_homology2=chr_index*2+1;

    int start_index1=homology_tad_start[chr_homology1]; 
    int start_index2=homology_tad_start[chr_homology2]; 

    Eigen::MatrixXf new_chromosome1=tads_center_chr[chr_homology1];
    Eigen::MatrixXf new_chromosome2=tads_center_chr[chr_homology2];
    
    homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+start_index1,homology_dis_matrix[thread_num].begin()+start_index2+tad_len); 

    chrom_move(new_chromosome1,move_index);
    float energy1=get_chrom_delta_E(thread_num,chr_homology1,new_chromosome1,tads_center_chr,homology_dis_matrix[thread_num]);
    float sphere_limition1=sphere_limit(new_chromosome1,alpha)-sphere_limit(tads_center_chr[chr_homology1],alpha);

    chrom_move(new_chromosome2,move_index);      
    float energy2=get_chrom_delta_E(thread_num,chr_homology2,new_chromosome2,tads_center_chr,homology_dis_matrix[thread_num]);
    float sphere_limition2=sphere_limit(new_chromosome2,alpha)-sphere_limit(tads_center_chr[chr_homology2],alpha);    
    //能量切换的时候e没有保留                                                                                                                                                       
    if(accept_move(energy1+sphere_limition1+energy2+sphere_limition2)){

        //update center      
        tads_center_chr[chr_homology1]=new_chromosome1;
        tads_center_chr[chr_homology2]=new_chromosome2;      

        //update bin_location
        chrom_move(connected_structure_chr[chr_homology1],move_index);
        chrom_move(connected_structure_chr[chr_homology2],move_index);

        // //update neighbor
        // update_chrom_dis_matrix(thread_num,chr_homology1); 
        // update_chrom_dis_matrix(thread_num,chr_homology2);

        // //update contacts
        // update_chrom_neighbor_contact(chr_homology1,thread_num,m);    
        // update_chrom_neighbor_contact(chr_homology2,thread_num,m);           
          
    }
 

    // string dis_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure/test.txt";
    // ofstream distance(dis_file);
    // for(int i=0;i<tad_len*2;i++){//2958
    //     for(int j=0;j<1500;j++){
    //         distance<<homology_dis_matrix_change[thread_num][i][j]<<'\t';
    //         // distance<<homology_dis_matrix[0][i][j]<<'\t';
    //     }
    //     distance<<'\n';
    // }
    // distance.close();        
}

void chromosome_rotate(int thread_num,vector<MatrixXi> &connected_structure_chr,vector<vector<MatrixXi>> &tads_chr,vector<MatrixXf> &tads_center_chr,vector<MatrixXi> &dir_match,int chr_homology,int m){

    Eigen::MatrixXf new_chromosome=tads_center_chr[chr_homology];
    uniform_int_distribution<int> unirotate(0,6); 

    int start_index=homology_tad_start[chr_homology]; 
    int tad_len=chr_tad_len[chr_homology/2];

    int rotate_index=unirotate(gen);
    chrom_rotate(new_chromosome,rotate_index);

    homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+start_index,homology_dis_matrix[thread_num].begin()+start_index+tad_len);          
    
    float energy=get_chrom_delta_E(thread_num,chr_homology,new_chromosome,tads_center_chr,homology_dis_matrix[thread_num]);
    float sphere_limition=sphere_limit(new_chromosome,alpha)-sphere_limit(tads_center_chr[chr_homology],alpha);
    
    if(accept_move(energy+sphere_limition)){

        //update center
        tads_center_chr[chr_homology]=new_chromosome;     

        //update bin_location
        chrom_rotate(connected_structure_chr[chr_homology],rotate_index);
        
        //update neighbor
        update_chrom_dis_matrix(thread_num,chr_homology); 

        //update contacts
        update_chrom_neighbor_contact(chr_homology,thread_num,m);            
    }

}

void chrom_rotate_move(int thread_num,vector<MatrixXi> &connected_structure_chr,vector<vector<MatrixXi>> &tads_chr,vector<MatrixXf> &tads_center_chr,vector<MatrixXi> &dir_match,int chr_homology,int tad_index,int m){

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
    MatrixXi old_tads;
    MatrixXi tads;
    int dir=dir_match[chr_homology](0,tad_index);

    uniform_int_distribution<int> uniconnect(0,7);
    int dir_change=uniconnect(gen); 
    while(dir==dir_change){
        dir_change=uniconnect(gen); 
    }

    //update m
    thread_m[thread_num]=m;     

    if(tad_index>=tad_mid_len ){

        //intial dis_matrix_change

        int start_index=homology_tad_start[chr_homology];
        homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+start_index+tad_index,homology_dis_matrix[thread_num].begin()+start_index+tad_len);        

        first_point=connected_structure_chr[chr_homology](start,all);
        tads_center_change=tads_center_chr[chr_homology].bottomRows(tad_len-tad_index);

        MatrixXf move_point(1,3);
        move_point<<float(first_point(0,0)),float(first_point(0,1)),float(first_point(0,2));  
        rotate_center(move_point,tads_center_change,dir_change);     

        //update homology_dis_matrix_change and caculate energy    
        float energy=get_delta_E(thread_num,chr_homology,tad_index,tads_center_change,tads_center_chr,homology_dis_matrix[thread_num]);

        //get sphere_limition
        float sphere_limition=sphere_limit(tads_center_change,alpha)-sphere_limit(tads_center_chr[chr_homology].bottomRows(tad_len-tad_index),alpha);            
 
        // if(energy!=0){cout<<1<<'\n';}
        //energy更新有问题吧
        if(accept_move(energy + sphere_limition)){

            //update coordinate 
            change_structure=connected_structure_chr[chr_homology].bottomRows(chr_bin_len-start);          
            rotate_center(first_point,change_structure,dir_change);  
            connected_structure_chr[chr_homology].bottomRows(chr_bin_len-start)=change_structure;
                                  
            //update tad_center   
            tads_center_chr[chr_homology].bottomRows(tad_len-tad_index)=tads_center_change;         

            //update dir_match  
            dir_match[chr_homology](0,tad_index)=dir_change;    

            //update neighbor
            update_dis_matrix(thread_num,chr_homology,tad_index,1);    

            //update contacts
            update_neighbor_contact(chr_homology,thread_num,m);

        }               
    }    
    else{

        //intial dis_matrix_change

        int start_index=homology_tad_start[chr_homology];
        homology_dis_matrix_change[thread_num].assign(homology_dis_matrix[thread_num].begin()+start_index,homology_dis_matrix[thread_num].begin()+start_index+tad_index+1);     

        first_point=connected_structure_chr[chr_homology](end,all);
        tads_center_change=tads_center_chr[chr_homology].topRows(tad_index+1);      

        MatrixXf move_point(1,3);
        move_point<<float(first_point(0,0)),float(first_point(0,1)),float(first_point(0,2));  
        rotate_center(move_point,tads_center_change,dir_change);          

        float energy=get_delta_E(thread_num,chr_homology,tad_index,tads_center_change,tads_center_chr,homology_dis_matrix[thread_num]);

        //get sphere_limition
        float sphere_limition=sphere_limit(tads_center_change,alpha)-sphere_limit(tads_center_chr[chr_homology].topRows(tad_index+1),alpha);            
        
        if(accept_move(energy + sphere_limition)){    

            //update coordinate 
            change_structure=connected_structure_chr[chr_homology].topRows(end);          
            rotate_center(first_point,change_structure,dir_change);  
            connected_structure_chr[chr_homology].topRows(end)=change_structure; 

            //update tad_center   
            tads_center_chr[chr_homology].topRows(tad_index+1)=tads_center_change;         

            //update dir_match  
            dir_match[chr_homology](0,tad_index)=dir_change;    

            //update neighbor
            update_dis_matrix(thread_num,chr_homology,tad_index,0); 

            //update contacts
            update_neighbor_contact(chr_homology,thread_num,m);              
                     
        }
    }

    //unlock
    mutexs[thread_num].unlock();    

    
}



#endif //Chrom_MOVES_H