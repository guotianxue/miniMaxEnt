#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
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


string outfiledir="";
string outfile_access="";
string corr_file="";
string out_polymer_filedir="";
string chr="";

bool initial_by_file=false;
bool data_by_input=false;
bool out_contact_only=false;
bool rough_refer=false;

int start_bin;
int end_bin;
int replicate_iter=-1;
int mark=0;
int reconstruct_chrom_mark=0;

//-------------外部参数-----------------

// float learning_rate=0.05;
float learning_rate=0.05;
int radius=20;//the radius of cell nucleus
int number_of_threads =48;//48
int update_steps =15;
long int burn_in_time = 1000000;
long int mc_moves_start = 2000;
float alpha=0;
chrono::duration<float> time_start;

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

float corr(vector<vector<float>> &combine){
    float sum1=0,sum2=0,sum3=0,mean1,mean2,value1,value2;
    int length=combine.size();
    for(int i=0;i<length;i++){
        for(int j=i+1;j<length;j++){
            sum1+=combine[i][j];
            sum2+=combine[j][i];
        }
    }

    mean1=2*sum1/length/length;
    mean2=2*sum2/length/length;

    sum1=0;
    sum2=0;

    for(int i=0;i<length;i++){
        for(int j=i+1;j<length;j++){
            value1=combine[i][j]-mean1;
            value2=combine[j][i]-mean2;

            sum3+=value1*value2*2;
            sum2+=value1*value1*2;
            sum1+=value2*value2*2;
        }
    }
    float r=sum3/sqrt(sum1*sum2);
    return r;

}

void homology_merge(int thread_num,int chr_index1,int chr_index2){

    int chr_len1,chr_len2,homology_start1,homology_start2,homology_start11,homology_start22,haploid_start1,haploid_start2;
    float contact1,contact2,contact3,contact4;

    chr_len1=chr_tad_len[chr_index1];
    chr_len2=chr_tad_len[chr_index2];
    homology_start1=homology_tad_start[chr_index1*2];
    homology_start11=homology_tad_start[chr_index1*2+1];
    homology_start2=homology_tad_start[chr_index2*2];
    homology_start22=homology_tad_start[chr_index2*2+1];    
    haploid_start1=haploid_tad_start[chr_index1];
    haploid_start2=haploid_tad_start[chr_index2];

    for(int i=0;i<chr_len1;i++){
        if(chr_index1==chr_index2){
            for(int j=i+1;j<chr_len2;j++){
                contact1=homology_contact_matrix[thread_num][homology_start1+i][homology_start2+j];
                contact2=homology_contact_matrix[thread_num][homology_start11+i][homology_start22+j];

                contact3=homology_contact_matrix[thread_num][homology_start1+i][homology_start22+j];
                contact4=homology_contact_matrix[thread_num][homology_start1+j][homology_start22+i];                
            
                haploid_contact_matrix[haploid_start1+i][haploid_start2+j]+=contact1+contact2+contact3+contact4; 
          
            }                        
        }
        else{
            for(int j=0;j<chr_len2;j++){
                contact1=homology_contact_matrix[thread_num][homology_start1+i][homology_start2+j];
                contact2=homology_contact_matrix[thread_num][homology_start11+i][homology_start2+j];
                contact3=homology_contact_matrix[thread_num][homology_start1+i][homology_start22+j];
                contact4=homology_contact_matrix[thread_num][homology_start11+i][homology_start22+j];                
            
                haploid_contact_matrix[haploid_start1+i][haploid_start2+j]+=contact1+contact2+contact3+contact4;                 
            }
        }
    }
}

void get_final_structure(int thread_num,int chr_index,MatrixXi &M,string out_polymer_filedir){
    int index;int multiply;vector <int> coordinate;
    int margin_abs;
    MatrixXi margin;
    // ofstream fout("/data/home/txguo/data_use/maxEnt/whole_chr/rough_big_structure/"+chr+"/"+to_string(thread_num)+".txt");

    // string outfile=outfiledir+"/polymer/chr"+chr+"_rough_"+to_string(thread_num)+".txt";

    ofstream fout(out_polymer_filedir+"/chr"+to_string(chr_index+1)+"_"+to_string(thread_num)+".txt");
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

void reconstruct_chrom(int chr_homology){

    //initial the dis to the center
    for(int thread_num=0;thread_num<number_of_threads;thread_num++){
        initial_tads_dis(thread_num,chr_homology);  
    }    

    //initial start time 
    auto start = chrono::high_resolution_clock::now();

    int chrom_tads_num=chr_tad_len[chr_homology/2];
           
    //cut genome level 2 chrom level
    cut_tads_refer(chr_homology);
    normalize_tads_reference(reference_tads_contact_s);       
    
    cout<<"============"<<'\n';
    cout<<"chr:"<<to_string(chr_homology)<<'\n';
    cout<<"m:"<<mc_moves_start<<'\n';
    cout<<"radius:"<<radius<<'\n';
    int multiple=-3;   //-3
    vector<vector<float>>(chrom_tads_num,vector<float>(chrom_tads_num,0)).swap(chr_E); 

    for (int n = 0; n<update_steps;n++) {   
        cout<<"============"<<'\n';
        //update alpha

        alpha=min(exp(n)*pow(10,multiple),0.5);
        // alpha=0.5;

        cout<<"alpha:"<<alpha<<'\n';
        cout<<"iter"<<n<<'\n';
        
        //initial
        for (auto l = 0; l < number_of_threads; l++) {    
            tads_center_m[l]= Eigen::MatrixXf::Ones(chrom_tads_num,chrom_tads_num)*(-1);

            //setzero homology_m
            vector<int> h_m(chr_homology_num,0);
            homology_m[l].swap(h_m);

            //setzero contact matrix
            vector<vector<float>> contact_matrix(chrom_tads_num,vector<float>(chrom_tads_num,0));
            homology_contact_matrix[l].swap(contact_matrix);   
 

        }
        //initial haploid
        // int haploid_tads_num=tads_num/2;   
        int haploid_tads_num=chrom_tads_num;
        vector<vector<float>>(haploid_tads_num,vector<float>(haploid_tads_num,0)).swap(haploid_contact_matrix);
        vector<int>(number_of_threads,0).swap(thread_m);       

        //mc moves

        cout<<"tads_move"<<'\n';
        mark=0;
        for (auto l = 0; l < number_of_threads; l++) {         
            threads[l] = thread(run_tads2chrom,l,mc_moves_start,chr_homology);//创建线程
        }
        for (auto &&l : threads) {
            l.join();
        }                           
    
        auto time1 = chrono::high_resolution_clock::now();
        chrono::duration<float> elapsed1 = time1 - start;
        cout<<"run_time:"<<elapsed1.count() << " seconds\n";  

        //merge 
        for(int thread_num=0;thread_num<number_of_threads;thread_num++){
            for(int i=0;i<chrom_tads_num;i++){
                for(int j=i+1;j<chrom_tads_num;j++){
                    haploid_contact_matrix[i][j]+=homology_contact_matrix[thread_num][i][j];
                }
            }
        }

        // //test
        // vector<vector<float>> contact_matrix(chrom_tads_num,vector<float>(chrom_tads_num,0));
        // homology_contact_matrix[0].swap(contact_matrix); 
        // for(int thread_num=0;thread_num<number_of_threads;thread_num++){
        //     for(int i=0;i<chrom_tads_num;i++){
        //         for(int j=i+1;j<chrom_tads_num;j++){
        //             homology_contact_matrix[0][i][j]+=1/pow(homology_dis_matrix[thread_num][i][j],2);
        //         }
        // }

        // string c_file="/data/home/txguo/code_final/reconstruct_chromosome_structure/contact.txt";
        // ofstream con(c_file);
        // for(int i=0;i<chrom_tads_num;i++){//2958
        //     for(int j=0;j<chrom_tads_num;j++){
        //         con<<homology_contact_matrix[0][i][j]<<'\t';
        //         // distance<<homology_dis_matrix[0][i][j]<<'\t';
        //     }
        //     con<<'\n';
        // }
        // con.close();               

        normalize_tads();   
        // //update energy & combine refer and real contact
        
        update_tads_energies(); 
        

        string e_file="/data/home/txguo/code_final/reconstruct_chromosome_structure/test.txt";
        ofstream energ(e_file);
        int chr_length=chrom_tads_num;
        for(int i=0;i<chr_length;i++){//2958
            for(int j=0;j<chr_length;j++){
                energ<<chr_E[i][j]<<'\t';
                // distance<<homology_dis_matrix[0][i][j]<<'\t';
            }
            energ<<'\n';
        }
        energ.close();                                        
 
        //输出homology_contact_matrix看一下结果
        // string dis_file="/data/home/txguo/code_final/reconstruct_chromosome_structure/dis.txt";
        // string contact_file="/data/home/txguo/data_use/maxEnt/whole_chr/chr_tads_structure/"+to_string(chr_homology+1)+"/chr"+to_string(chr_homology+1)+"_contact.txt";
        string contact_file=outfiledir+"/"+to_string(chr_homology+1)+"/chr"+to_string(chr_homology+1)+"_contact.txt";
        ofstream distance(contact_file);
        for(int i=0;i<chrom_tads_num;i++){//2958
            for(int j=0;j<chrom_tads_num;j++){
                distance<<haploid_contact_matrix[i][j]<<'\t';
                // distance<<homology_dis_matrix_change[0][i][j]<<'\t';
            }
            distance<<'\n';
        }
        distance.close();                           

        // get corr
        float cc=corr(haploid_contact_matrix);
        cout<<"corr:"<<cc<<'\n';                        

        //output center
        for(int i=0;i<number_of_threads;i++){
            // string p_file="/data/home/txguo/code_final/reconstruct_chromosome_structure/polymer"+to_string(i)+".txt";
            // string c_file="/data/home/txguo/data_use/maxEnt/whole_chr/chr_tads_structure/"+to_string(chr_homology+1)+"/center"+to_string(i)+".txt";
            string c_file=outfiledir+"/"+to_string(chr_homology+1)+"/center"+to_string(i)+".txt";
            ofstream Center(c_file);

            int length=tads_center_chr[i][chr_homology].rows();
            MatrixXf &tmp = tads_center_chr[i][chr_homology];

            for (int l=0;l<length;l++){
                Center<<chr_homology<<'\t'<<tmp(l,0)<<'\t'<<tmp(l,1)<<'\t'<<tmp(l,2)<<'\n';
            }

            Center.close();              
        }

                 

        // // get corr
        // float c=corr(haploid_contact_matrix);
        // cout<<"inter_corr"<<'\n';
        // cout<<"corr:"<<c<<'\n';  

        //output polymer for plot
        if(n%5==0){
            for(int i=0;i<5;i++){
                // string p_file="/data/home/txguo/data_use/maxEnt/whole_chr/chr_tads_structure/"+to_string(chr_homology+1)+"/iter/chr"+to_string(chr_homology+1)+"_iter"+to_string(n)+"_"+to_string(i)+".txt";
                string p_file=outfiledir+"/"+to_string(chr_homology+1)+"/iter/chr"+to_string(chr_homology+1)+"_iter"+to_string(n)+"_"+to_string(i)+".txt";
                ofstream P(p_file);

                int length=chr_len[chr_homology/2];
                MatrixXi &tmp = connected_structure_chr[i][chr_homology];
                
                for (int l=0;l<length;l++){
                    P<<chr_homology<<'\t'<<tmp(l,0)<<'\t'<<tmp(l,1)<<'\t'<<tmp(l,2)<<'\n';
                }

                P.close();     

            }
        }


        //output single
        string file1;
        for(int k=0;k<8;k++){
            // file1="/data/home/txguo/data_use/maxEnt/whole_chr/chr_tads_structure/"+to_string(chr_homology+1)+"/homology/chr"+to_string(chr_homology+1)+"_"+to_string(k)+"_contact.txt";
            file1=outfiledir+"/"+to_string(chr_homology+1)+"/single/chr"+to_string(chr_homology+1)+"_"+to_string(k)+"_contact.txt";
            ofstream contact1(file1);
            for(int i=0;i<chr_length;i++){//2958
                for(int j=0;j<chr_length;j++){
                    contact1<<homology_contact_matrix[k][i][j]<<'\t';
                    // distance<<homology_dis_matrix[0][i][j]<<'\t';
                }
                contact1<<'\n';
            }       
            contact1.close();        
        }

           
    }   

    //output polymer
    for(int i=0;i<number_of_threads;i++){
        // string p_file="/data/home/txguo/data_use/maxEnt/whole_chr/chr_tads_structure/"+to_string(chr_homology+1)+"/polymer"+to_string(i)+".txt";  
        string p_file=outfiledir+"/"+to_string(chr_homology+1)+"/polymer"+to_string(i)+".txt";
        ofstream P(p_file);

        int length=chr_len[chr_homology/2];
        MatrixXi &tmp = connected_structure_chr[i][chr_homology];
        
        for (int l=0;l<length;l++){
            P<<chr_homology<<'\t'<<tmp(l,0)<<'\t'<<tmp(l,1)<<'\t'<<tmp(l,2)<<'\n';
        }

        P.close();     

        // string out_polymer_filedir="/data/home/txguo/data_use/maxEnt/whole_chr/chr_tads_structure/"+to_string(chr_homology+1)+"/polymer";
        string out_polymer_filedir=outfiledir+"/"+to_string(chr_homology+1)+"/polymer";
        get_final_structure(i,chr_homology,connected_structure_chr[i][chr_homology],out_polymer_filedir);         
    }            

    // //output polymer
    // for(auto l=0;l<number_of_threads;l++){
    //     get_final_structure(l,out_polymer_filedir);
    // }

}

void reconstruct(){
    //---------initial mutexs------------

    // string tad_file="/data/home/txguo/data_use/maxEnt/whole_chr/divide_tad/chr"+chr+"_tad_divide.txt";   
    get_tadset(tads_bound);

    // initial_chrom_point();
    initial_chrom_point_zero();

    string file="/data/home/txguo/data_use/maxEnt/whole_chr/small_tads/polymer"; 
    
    // initial tads file
    for (auto l = 0; l < number_of_threads; ++l) {
        threads[l] = thread(initial_tads_byfile,l,ref(tads_chr[l]),ref(connected_structure_chr[l]),file);
    }
    for (auto &&l : threads) {
        l.join();
    }

    initial_tads_data();  

    //normalize reference
    initial_tad_contact_refer() ;

    cout<<"initial_refer_over"<<'\n';

    // for(int chr_homology=0;chr_homology<1;chr_homology++){

    int chr_homology=stoi(chr);
    reconstruct_chrom_mark=1;
    tads_num=chr_tad_len[chr_homology/2];  
    // radius=pow(tads_num/2,0.5);
    // radius=ceil(pow(tads_num*3,0.33));
    // radius=20;
    // cout<<"radius:"<<radius<<'\n';
    reconstruct_chrom(chr_homology);
    // }      

    
}