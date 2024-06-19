#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <random>
#include <string>
#include <thread>
#include <iomanip>
#include <chrono>
#include "initial_polymer.h"
#include "Moves_liner.h"
#include "Inverse_Monte_Carlo.h"

using namespace std;
using namespace Eigen;


const int refer_length=12463;//chr1
long int mc_moves;

string chr="1";
// string outfile="/data/home/txguo/data_use/maxEnt/gm12878_tad/contacts";
string outfile="";
string inputfiledir="";

bool initial_tads=false;
bool data_by_input=false;
bool out_contact_only=false;
bool rough_refer=false;

int start_bin;
int end_bin;
int replicate_iter=-1;
int beta=0;
//-------------外部参数-----------------

float learning_rate=0.1;//0.01;
int radius=55;//the radius of cell nucleus
int number_of_threads = 48;//48
// int update_steps =20;
int update_steps =500;
long int burn_in_time = 100000;//10000;
long int mc_moves_start = 1;
// vector<int> chr_r_list={15,14,14,14,13,13,13,12,12,12,11,11,11,11,11,9,9,9,7,7,6,6};
vector<int> chr_r_list={14,14,14,14,13,13,13,12,12,12,11,11,11,11,11,9,9,9,7,7,6,6};
// 12,10,10,10,9,9,9,9,9,7,7,7,6,6,6,6};

int mask_bins=2;

//neigbor 程序参数
int threshold=200*6;//if bin dis >=200 ,deal with neighbor in 3d
int neighbor_dis=5;//if dis>neighbor_dis,regard as a neighbor contact

float neighbor_energy_weight=0.0001; 
float neighbor_contact_weight=0.3;
//-------------------------------------

uniform_int_distribution<int> unisite;//(0,pol_length-1);
uniform_int_distribution<int> unipath(0,5);

uniform_real_distribution<float> unif(0.0,1);//                            产生均匀分布在区间 [a, b) 上的随机浮点值 i  P=1/(b − a)
uniform_int_distribution<int> unimove(0,3);// 3选择四种移动中的一种              生成随机整数值 i ，均匀分布于闭区间 [a, b] P=1/(b − a+1)
uniform_int_distribution<int> leftOright(0,1);
uniform_int_distribution<int> uniradius(-1*radius,radius);

vector<thread> threads(number_of_threads);

void structure_initial(){

    //---------------------- initial ----------------------------
    //initial

    for (auto l = 0; l < number_of_threads; l++) {    
        chrom_total_contacts[l]=Eigen::MatrixXd::Zero(chr_num,chr_num);
    }    
    chrom_final_contacts=Eigen::MatrixXd::Zero(chr_num,chr_num);

    initial_chrom_energy();
    initial_chrom_polymer(); 
    initial_chrom_refer();
    normalize_reference(reference_chrom_contact);

    return ;
}

//------------ get corr ---------
double mean(MatrixXd Matrix){
    return Matrix.sum()/(chr_num*chr_num);
}

float corr2(MatrixXd Matrix1,MatrixXd Matrix2){

    double mean1=mean(Matrix1);double mean2=mean(Matrix2);
    MatrixXd mean1_mat=MatrixXd::Constant(chr_num,chr_num, mean1);
    MatrixXd mean2_mat=MatrixXd::Constant(chr_num,chr_num, mean2);

    Matrix1=Matrix1-mean1_mat;
    Matrix2=Matrix2-mean2_mat;    
    
    float r=(Matrix1.cwiseProduct(Matrix2)).sum()/sqrt((Matrix1.cwiseProduct(Matrix1)).sum() * (Matrix2.cwiseProduct(Matrix2)).sum());
    return r;
}

void reconstruct(){

    // show parameter

    cout<<"#######################"<<'\n';
    cout<<"parameter:"<<'\n';
    cout<<"chr:"<<chr<<'\n';
    cout<<"radius:"<<'\t'<<radius<<'\n';
    cout<<"number_of_threads:"<<'\t'<<number_of_threads<<'\n';
    cout<<"update_steps:"<<'\t'<<update_steps<<'\n';
    cout<<"burn_in_time:"<<'\t'<<burn_in_time<<'\n';

    structure_initial();

    auto start = chrono::high_resolution_clock::now();

    for (int n =0; n<update_steps;n++) { //do iterative update scheme

        //------------------- initial ---------------------
        cout<<"---------------"<<'\n';
        reduce_count=0,add_count=0,kink_count=0,loop_count=0,skip_count=0;

        cout<<"iter_step:"<< n<<"\n";
        mc_moves = mc_moves_start*sqrt(n+10)/sqrt(10); //number of MC moves grows with each iteration (implicitly converted to long int)

        //reset contact frequencies before starting new forward round
        chrom_final_contacts=Eigen::MatrixXd::Zero(chr_num,chr_num);        
        for (auto l = 0; l < number_of_threads; l++) {    
            chrom_total_contacts[l]=Eigen::MatrixXd::Zero(chr_num,chr_num);
        }  

        //run forward simulation
        for (auto l = 0; l < number_of_threads; l++) {
            threads[l] = thread(run, mc_moves,l);
        }
        for (auto &&l : threads) {
            l.join();
        }
    
        //add up contacts from threads
        for (int i = 0; i < number_of_threads; i++) {
            chrom_final_contacts += chrom_total_contacts[i];
        }   

        //normalize contact frequencies
        normalize();           

        // //update energies
        update_energy_chrom(reference_chrom_contact);   

        beta++;                       

        auto time1 = chrono::high_resolution_clock::now();
        chrono::duration<float> elapsed1 = time1 - start;
        cout<<"run_time:"<<elapsed1.count() << " seconds\n";

        float corr=corr2(chrom_final_contacts+chrom_final_contacts.transpose(),reference_chrom_contact+reference_chrom_contact.transpose());
        cout<<"corr:"<<corr<<'\n';           

        std::ofstream out1("./test.txt");
        out1<<reference_chrom_contact+chrom_final_contacts.transpose();   
        out1.close();       

        std::ofstream out2("./e.txt");
        out2<<Interaction_E_chrom;
        out2.close();  

         
        for(int thread=0;thread<48;thread++) {
            string file="/data/home/txguo/data_use/maxEnt/whole_chr/chrom_initial/chrom_initial_"+to_string(thread)+".txt";
            int m=thread;
            std::ofstream polymer_res(file); 
            for(int k=0;k<polymer_chrom[m].size();k++){
                polymer_res<<polymer_chrom[m][k][0]<<'\t'<<polymer_chrom[m][k][1]<<'\t'<<polymer_chrom[m][k][2]<<'\n';  
            }   
            polymer_res.close();          
        }

                                 
    }                         
}



