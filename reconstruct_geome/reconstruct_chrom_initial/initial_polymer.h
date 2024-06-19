#ifndef Inverse_ploymer_Initialize_h
#define Inverse_ploymer_Initialize_h


#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <string>
#include "global.h"
#include "Energy_changes_liner.h"
using namespace std;
using namespace Eigen;


int tad_bin_num;
int pol_length;//monomers的数量
int chr_len;
int chr_num=22;

static mt19937_64 gen(time(0));//mt19937是伪随机数产生器，用于产生高性能的随机数

vector<vector<vector<int>>> polymer_chrom(number_of_threads);

//energy initial 
Eigen::MatrixXd Interaction_E_chrom;

//initial reference
Eigen::MatrixXd reference_chrom_contact;

//initial contact 
vector<Eigen::MatrixXd> chrom_total_contacts(number_of_threads);
Eigen::MatrixXd chrom_final_contacts;


//------------------------------polymer初始化------------------------------------------------
bool check_nucleus(int x,int y,int z){
    return (pow(x,2)+pow(y,2)+pow(z,2)<=pow(radius,2));
}

void initial_chrom_polymer(){
    
    // uniform_int_distribution<int> cordinate_choose(-15,15);
    uniform_int_distribution<int> cordinate_choose(-30,30);    
    int x,y,z;
    for (int k=0;k<number_of_threads;k++){
        for(int i=0;i<chr_num*2;i++){
            x=cordinate_choose(gen);
            y=cordinate_choose(gen);
            z=cordinate_choose(gen);
            polymer_chrom[k].push_back({x,y,z});
        }        
    }


}

//-----------------------tad,refer contact 和tad contact处理--------------------------

void initial_chrom_refer(){

    // initial the refer data

    string refer_file="/data/home/txguo/data_use/maxEnt/reference/reference_use/genome_contact.mtx";

    int sitei;
    int sitej;
    float contact_num;

    SparseMatrix<double> refer ;
    Eigen::loadMarket(refer,refer_file);
    reference_chrom_contact=MatrixXd(refer);
}


void initial_chrom_energy(){

    Interaction_E_chrom=Eigen::MatrixXd::Zero(chr_num,chr_num);

}


#endif