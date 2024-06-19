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
#include "mutex_use.h"
#include "initial_polymer.h"
#include "Moves_liner.h"
#include "Inverse_Monte_Carlo.h"

// g++ -std=c++11 -o test -I /data/software/HiC/vscode_package/eigen-3.4.0 -I /data/software/HiC/vscode_package/boost_1_79_0 ./test.cpp -pthread

using namespace std;
using namespace Eigen;


const int refer_length=12463;//chr1
long int mc_moves;

string chr="2";
string outfile="/data/home/txguo/data_use/maxEnt/gm12878_tad/contacts";

bool initial_by_file=false;
bool data_by_input=false;
bool out_contact_only=false;
bool rough_refer=false;

int start_bin;
int end_bin;
int replicate_iter=-1;

//-------------外部参数-----------------

float learning_rate=0.15;
int radius=20;//the radius of cell nucleus
int number_of_threads = 48;//48
int update_steps =5;
long int burn_in_time = 1000000;
long int mc_moves_start = 3000;
float alpha=0;

//neigbor 程序参数
int threshold=200*6;//if bin dis >=200 ,deal with neighbor in 3d
int neighbor_dis=5;//if dis>neighbor_dis,regard as a neighbor contact
// float neighbor_energy_weight=0.0001; 
float neighbor_energy_weight=0.0001; 
float neighbor_contact_weight=0.3;
//-------------------------------------

uniform_int_distribution<int> unisite;//(0,pol_length-1);
uniform_int_distribution<int> unipath(0,6);
uniform_int_distribution<int> uniTranslate(0,17);
uniform_int_distribution<int> uniThread(0,number_of_threads-1);

uniform_real_distribution<float> unif(0.0,1);//                            产生均匀分布在区间 [a, b) 上的随机浮点值 i  P=1/(b − a)
// uniform_int_distribution<int> unimove(0,2);// 选择三种移动中的一种              生成随机整数值 i ，均匀分布于闭区间 [a, b] P=1/(b − a+1)
uniform_int_distribution<int> unimove(0,9);
uniform_int_distribution<int> leftOright(0,1);
uniform_int_distribution<int> uniradius(-1*radius,radius);

vector<thread> threads(number_of_threads);

//--------------------
float mean(MatrixXf Matrix){
    int length=Matrix.rows();
    return Matrix.sum()/(length*length);
}

float corr(MatrixXf Matrix1,MatrixXf Matrix2){
    int length=Matrix1.rows();
    double mean1=mean(Matrix1);double mean2=mean(Matrix2);
    MatrixXf mean1_mat=MatrixXf::Constant(length,length, mean1);
    MatrixXf mean2_mat=MatrixXf::Constant(length,length, mean2);

    Matrix1=Matrix1-mean1_mat;
    Matrix2=Matrix2-mean2_mat;    
    
    float r=((Matrix1.cwiseProduct(Matrix2)).sum())/sqrt((Matrix1.cwiseProduct(Matrix1)).sum() * (Matrix2.cwiseProduct(Matrix2)).sum());
    return r;
}

void get_final_structure(int thread_num){
    connected_structure[thread_num];
    int index;int multiply;vector <int> coordinate;
    int margin_abs;
    MatrixXi margin;
    MatrixXi &M=connected_structure[thread_num];
    // ofstream fout("/data/home/txguo/data_use/maxEnt/whole_chr/rough_big_structure/"+chr+"/"+to_string(thread_num)+".txt");
        ofstream fout("/data/home/txguo/data_use/maxEnt/whole_chr/rough_big_structure/short/"+to_string(thread_num)+".txt");
    for (int i=0;i<M.rows()-1;++i){
        index=i*6;
        for(int j=0;j<M.cols();++j){
            fout << M(i,j) <<"\t";
        }
        fout<<index<<'\n';

        margin=M.row(i+1)-M.row(i);
        margin_abs=(margin.array().abs()).sum();
        coordinate={M(i,0),M(i,1),M(i,2)};

        if(margin_abs==0){
            fout << coordinate[0]-1 <<"\t"<< coordinate[1] <<"\t"<< coordinate[2]<<"\t"<< index+3 <<"\n";
        }
        else if (margin_abs==2){
            for (int k=0;k<3;k++){
                if (margin(0,k)!=0 ){
                    coordinate[k]+=margin(0,k)/abs(margin(0,k));
                    fout << coordinate[0] <<"\t"<< coordinate[1] <<"\t"<< coordinate[2]<<"\t"<< index+3 <<"\n";         
                    break;                
                }       
            }
        }
        else{
            for(int k=0;k<3;k++){
                while(margin(0,k)!=0 and (margin.array().abs()).sum()!=1){
                    index+=1;
                    multiply+=1;
                    for(int num=0;num<3;num++){
                        if(num==k){
                            coordinate[num]=coordinate[num]+margin(0,k)/abs(margin(0,k)) ;
                            fout<<coordinate[num]<<"\t";
                            }
                        else{fout << coordinate[num] <<"\t";}
                    }
                    fout << index <<"\n";
                    margin(0,k)-=margin(0,k)/abs(margin(0,k));
                    
                }
            }
        }              
    }
    //fill the last bin
    int last=M.rows()-1;
    index=last*6;
    fout << M(last,0) <<"\t"<< M(last,1) <<"\t"<< M(last,2) <<"\t"<< index <<"\n";
    fout << M(last,0)+1 <<"\t"<< M(last,1) <<"\t"<< M(last,2) <<"\t"<< index+3 <<"\n";
    fout << M(last,0)+2 <<"\t"<< M(last,1) <<"\t"<< M(last,2) <<"\t"<< index+4 <<"\n";
    fout << M(last,0)+3 <<"\t"<< M(last,1) <<"\t"<< M(last,2) <<"\t"<< index+5 <<"\n";      
    
    fout.close();      
}

int main(){
    //---------initial mutexs------------

    get_tadset("/data/home/txguo/tmp.txt",tads_bound);
    string file="/data/home/txguo/data_use/maxEnt/parameter_access/radius/25bins/r_5_n_0/polymer/2";

    cout<<"go"<<'\n';    

    // initial tads file
    for (auto l = 0; l < number_of_threads; ++l) {
        threads[l] = thread(initial_tads_byfile,l,ref(tads[l]),ref(connected_tads[l]),ref(connected_structure[l]),file,ref(tads_bound));
    }
    for (auto &&l : threads) {
        l.join();
    }


    initial_tads_data();

    // //normalize reference
    initial_rough_refer(chr,500) ;
 
    normalize_tads_reference();
    
    cout<<"go"<<'\n';
    int multiple=-3;   

    for (int n = 0; n<update_steps;n++) {   
        cout<<"============"<<'\n';
        //update alpha

        alpha=min(exp(n)*pow(10,multiple),1.00);
        cout<<"alpha:"<<alpha<<'\n';
        
        //initial
        for (auto l = 0; l < number_of_threads; l++) {    
            tads_center_m[l]= Eigen::MatrixXf::Ones(tads_bound.size(),tads_bound.size())*(-1);
            tads_total_contacts[l]=Eigen::MatrixXf::Zero(tads_bound.size(),tads_bound.size());
        }

        //mc moves
        for (auto l = 0; l < number_of_threads; l++) {         
            threads[l] = thread(run_tads,l,mc_moves_start);//创建线程
        }
        for (auto &&l : threads) {
            l.join();
        }
        
        // read in remaining contacts at the end of forward simulation
        for (auto l = 0; l < number_of_threads; l++) {
            for (auto elem : tads_contact[l]) { //add the contacts remaining at the end of the simulation
                tads_total_contacts[l](tad_bin_dict[elem.first.first],tad_bin_dict[elem.first.second]) += (mc_moves_start - elem.second);       
            }
            
            //clear old contact
            tads_contact[l].clear();   
            threads_m[l]=0;   
        } 

        //add up contacts from threads
        tads_final_contacts=Eigen::MatrixXf::Zero(tads_bound.size(),tads_bound.size());
        int length=tads_bound[tads_bound.size()-1][1]+1-tads_bound[0][0];     
        for (int i = 0; i < number_of_threads; i++) {
            tads_final_contacts += tads_total_contacts[i];
        }   

        normalize_tads();

        //update energy
        update_tads_energies();

        //get corr
        float c=corr(tads_final_contacts+tads_final_contacts.transpose(),reference_tads_contact+reference_tads_contact.transpose());
        
        cout<<"iter:"<<n<<'\n';
        cout<<"corr:"<<c<<'\n';  

        ofstream fout("/data/home/txguo/code/code_use_linux_pure_5.5/tmp.txt");
        MatrixXf M = tads_final_contacts+reference_tads_contact.transpose();
        for (int i=0;i<M.rows();++i){
            for(int j=0;j<M.rows();++j){
                fout << M(i,j) <<"\t";
            }
            fout<<'\n';
        }
        fout.close();  

        ofstream Cout("/data/home/txguo/code/code_use_linux_pure_5.5/test.txt");
        Cout << Interaction_E_tads;
        Cout.close();    

        ofstream Sout("/data/home/txguo/code/code_use_linux_pure_5.5/structure.txt");
        Sout << connected_structure[0];
        Sout.close();           

    }   
    for(auto l=0;l<number_of_threads;l++){
        get_final_structure(l);
    }

    // //output polymer structure
    // for (auto l = 0; l < number_of_threads; l++) {         
    //     threads[l] = thread(get_final_structure,l);//创建线程
    // }
    // for (auto &&l : threads) {
    //     l.join();
    // }            
}