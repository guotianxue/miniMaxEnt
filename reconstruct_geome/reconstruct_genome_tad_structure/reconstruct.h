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

string chr="2";
string outfiledir="";
string outfile_access="";
string corr_file="";
string out_polymer_filedir="";

bool initial_by_file=false;
bool data_by_input=false;
bool out_contact_only=false;
bool rough_refer=false;

int start_bin;
int end_bin;
int replicate_iter=-1;
int mark=0;
int structure_initial_mark=0;//mark=0 initial randomly;  mark=1 initial by file;

//-------------外部参数-----------------

// float learning_rate=0.05;
float learning_rate=0.5;
int radius=60;//45;//the radius of cell nucleus
int number_of_threads =48;//48
int update_steps =1;
long int burn_in_time = 1000000;
long int mc_moves_start =10;
float alpha=0;

//neigbor 程序参数g
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

void homology_merge_dis(int thread_num,int chr_index1,int chr_index2){

    int homology_1,homology_2,chr_len1,chr_len2,homology_start1,homology_start2,homology_start11,homology_start22,haploid_start1,haploid_start2;
    float dis1,dis2,dis3,dis4;

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

                dis1=homology_dis_matrix[thread_num][homology_start1+i][homology_start2+j];
                dis2=homology_dis_matrix[thread_num][homology_start11+i][homology_start22+j];

                dis3=homology_dis_matrix[thread_num][homology_start1+i][homology_start22+j];
                dis4=homology_dis_matrix[thread_num][homology_start1+j][homology_start22+i];
                   
                haploid_contact_matrix[haploid_start1+i][haploid_start2+j]+=1/dis1/dis1+1/dis2/dis2+1/dis3/dis3+1/dis4/dis4;
            }
        }
        else{
            for(int j=0;j<chr_len2;j++){
                dis1=homology_dis_matrix[thread_num][homology_start1+i][homology_start2+j];
                dis2=homology_dis_matrix[thread_num][homology_start11+i][homology_start2+j];
                dis3=homology_dis_matrix[thread_num][homology_start1+i][homology_start22+j];
                dis4=homology_dis_matrix[thread_num][homology_start11+i][homology_start22+j];                   
                haploid_contact_matrix[haploid_start1+i][haploid_start2+j]+=1/dis1/dis1+1/dis2/dis2+1/dis3/dis3+1/dis4/dis4;
            }           
        }
    }
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

    // cout<<"test:"<<homology_contact_matrix[0][0][1]<<'\n';

    for(int i=0;i<chr_len1;i++){
        if(chr_index1==chr_index2){
            for(int j=i+1;j<chr_len2;j++){
                contact1=homology_contact_matrix[thread_num][homology_start1+i][homology_start2+j];
                contact2=homology_contact_matrix[thread_num][homology_start11+i][homology_start22+j];

                contact3=homology_contact_matrix[thread_num][homology_start1+i][homology_start22+j];
                contact4=homology_contact_matrix[thread_num][homology_start1+j][homology_start22+i];                
            
                haploid_contact_matrix[haploid_start1+i][haploid_start2+j]+=contact1+contact2+contact3+contact4; 

                // if(contact1<0 || contact2<0 || contact3<0 || contact4<0) {
                //     cout<<i<<'\t'<<j<<'\t'<<contact1<<'\t'<<contact2<<'\t'<<contact3<<'\t'<<contact4<<'\n';
                //     cout<<haploid_contact_matrix[haploid_start1+i][haploid_start2+j]<<'\n';
                // }               
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

void get_final_structure(int thread_num,int chr_homology,string out_polymer_filedir){
    connected_structure[thread_num];
    int index;int multiply;vector <int> coordinate;
    int margin_abs;
    MatrixXi margin;
    MatrixXi &M=connected_structure[thread_num];
    // ofstream fout("/data/home/txguo/data_use/maxEnt/whole_chr/rough_big_structure/"+chr+"/"+to_string(thread_num)+".txt");

    // string outfile=outfiledir+"/polymer/chr"+chr+"_rough_"+to_string(thread_num)+".txt";
    ofstream fout(out_polymer_filedir+"/chr"+chr+"_"+to_string(thread_num)+".txt");
    for (int i=0;i<M.rows()-1;++i){
        index=i*6;
        fout<<chr_homology<<'\t';
        for(int j=0;j<M.cols();++j){
            fout << M(i,j) <<"\t";
        }
        fout<<index<<'\n';

        margin=M.row(i+1)-M.row(i);
        margin_abs=(margin.array().abs()).sum();
        coordinate={M(i,0),M(i,1),M(i,2)};

        if(margin_abs==0){
            fout<<chr_homology<<'\t'<< coordinate[0]-1 <<"\t"<< coordinate[1] <<"\t"<< coordinate[2]<<"\t"<< index+3 <<"\n";
        }
        else if (margin_abs==2){
            for (int k=0;k<3;k++){
                if (margin(0,k)!=0 ){
                    coordinate[k]+=margin(0,k)/abs(margin(0,k));
                    fout<<chr_homology<<'\t'<< coordinate[0] <<"\t"<< coordinate[1] <<"\t"<< coordinate[2]<<"\t"<< index+3 <<"\n";         
                    break;                
                }       
            }
        }
        else{
            for(int k=0;k<3;k++){
                while(margin(0,k)!=0 and (margin.array().abs()).sum()!=1){
                    index+=1;
                    multiply+=1;
                    fout<<chr_homology<<'\t';
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
    fout<<chr_homology<<'\t'<< M(last,0) <<"\t"<< M(last,1) <<"\t"<< M(last,2) <<"\t"<< index <<"\n";
    fout<<chr_homology<<'\t'<< M(last,0)+1 <<"\t"<< M(last,1) <<"\t"<< M(last,2) <<"\t"<< index+3 <<"\n";
    fout<<chr_homology<<'\t'<< M(last,0)+2 <<"\t"<< M(last,1) <<"\t"<< M(last,2) <<"\t"<< index+4 <<"\n";
    fout<<chr_homology<<'\t'<< M(last,0)+3 <<"\t"<< M(last,1) <<"\t"<< M(last,2) <<"\t"<< index+5 <<"\n";      
    
    fout.close();      
}

void reconstruct(){
    //---------initial mutexs------------

    // string tad_file="/data/home/txguo/data_use/maxEnt/whole_chr/divide_tad/chr"+chr+"_tad_divide.txt";   
    get_tadset(tads_bound);

    initial_chrom_point();

    string file="/data/home/txguo/data_use/maxEnt/whole_chr/small_tad/polymer"; 

    structure_initial_mark=1;

    // initial tads file
    for (auto l = 0; l < number_of_threads; ++l) {
        threads[l] = thread(initial_tads_byfile,l,ref(tads_chr[l]),ref(connected_structure_chr[l]),file);
    }
    for (auto &&l : threads) {
        l.join();
    }

    //normalize reference
    initial_tad_contact_refer() ;
 
    normalize_tads_reference();      

    initial_tads_data();
     
    
    cout<<"============"<<'\n';
    cout<<"chr:"<<chr<<'\n';
    cout<<"radius:"<<radius<<'\n';
    cout<<"m:"<<mc_moves_start<<'\n';
    int multiple=-3;   //-3

    auto start = chrono::high_resolution_clock::now();
    vector<vector<string>> corr_list;
    for (int n = 0; n<update_steps;n++) {   
        cout<<"============"<<'\n';
        //update alpha

        alpha=min(exp(n)*pow(10,multiple),0.5);
        // alpha=0;
        cout<<"alpha:"<<alpha<<'\n';
        cout<<"iter"<<n<<'\n';
        
        //initial
        for (auto l = 0; l < number_of_threads; l++) {    
            tads_center_m[l]= Eigen::MatrixXf::Ones(tads_bound.size(),tads_bound.size())*(-1);

            //setzero homology_m
            vector<int> h_m(chr_homology_num,0);
            homology_m[l].swap(h_m);

            //setzero contact matrix
            vector<vector<float>> contact_matrix(tads_num,vector<float>(tads_num,0));
            homology_contact_matrix[l].swap(contact_matrix);    
        }
        (vector<vector<float>>(tads_num/2,vector<float>(tads_num/2,0))).swap(haploid_contact_matrix);
        vector<int>(number_of_threads,0).swap(thread_m);

        //mc moves


        // if(n>=10){            
        //     cout<<"tads_move"<<'\n';
        //     if(n==10){
        //         initial_E_By_structure();
        //         (vector<vector<float>>(tads_num/2,vector<float>(tads_num/2,0))).swap(haploid_contact_matrix);
        //     }
        //     mark=0;
        //     for (auto l = 0; l < number_of_threads; l++) {         
        //         threads[l] = thread(run_tads,l,mc_moves_start);//创建线程
        //     }
        //     for (auto &&l : threads) {
        //         l.join();
        //     }                                                           
        // }
        // else{
        //     cout<<"chrom_move"<<'\n';
        //     mark=1;            
        //     for (auto l = 0; l < number_of_threads; l++) {         
        //         threads[l] = thread(run_chrom_tads,l,mc_moves_start);//创建线程
        //     }
        //     for (auto &&l : threads) {
        //         l.join();
        //     }                         
        // }        


        cout<<"tads_move"<<'\n';
        mark=0;
        for (auto l = 0; l < number_of_threads; l++) {         
            threads[l] = thread(run_tads,l,mc_moves_start);//创建线程
        }
        for (auto &&l : threads) {
            l.join();
        }       

        // cout<<"chrom_move"<<'\n';
        // mark=1;            
        // for (auto l = 0; l < number_of_threads; l++) {         
        //     threads[l] = thread(run_chrom_tads,l,mc_moves_start);//创建线程
        // }
        // for (auto &&l : threads) {
        //     l.join();
        // }                            

        
        auto time1 = chrono::high_resolution_clock::now();
        chrono::duration<float> elapsed1 = time1 - start;
        cout<<"run_time:"<<elapsed1.count() << " seconds\n";  

        //merge homology contact matrix 2 haploid And add up contacts from threads

        // for(int k=0;k<tads_num;k++){
        //     for(int j=0;j<tads_num;j++){

        //         if(homology_contact_matrix[0][k][j]<0){
        //             cout<<k<<'\t'<<j<<'\t'<<homology_contact_matrix[0][k][j]<<'\n';
        //         }
        //     }
        // }           


        for (int l = 0; l < number_of_threads; l++)
        {
            for(int i=0;i<chr_num;i++){
                for(int j=i;j<chr_num;j++){                   
                    homology_merge(l,i,j);  
                }     
            }     
        }
     
        // if(n>10 && n%2!=1){
        //     normalize_tads();   
        //     //update energy & combine refer and real contact
        //     update_tads_energies();                            
              
        // }
        // else if(n<10){
        //     normalize_tads();     
        //     update_chrom_energies2();                     
        // }
        // else{
        //     normalize_chrom(); 
        //     // normalize_tads(); 
        //     // //update energy & combine refer and real contact
        //     // update_chrom_energies();   
        //     update_chrom_energies();                                    
        // }        

        // if(n<10){
        //     cout<<"chrom"<<'\n';
        //     normalize_tads(); 
        //     // // //update energy & combine refer and real contact
        //     update_chrom_energies2();             
        // }
        // else{
        //     cout<<"tad"<<'\n';
        //     // normalize_chrom(); 
        //     normalize_tads(); 
        //     // // //update energy & combine refer and real contact
        //     // update_chrom_energies();    
        //     update_tads_energies(); 
        // }

        // normalize_chrom(); 
        // update_chrom_energies();   
        normalize_tads();
        update_tads_energies();    
        // normalize_chrom();
        // update_chrom_energies();         

        // //test
        // int thread_num=0;
        // int homology=1;
        // int homology2=3;
        // int chr_len=chr_tad_len[homology/2];
        // int homology_start=homology_tad_start[homology];
        // int homology_end=homology_start+chr_len;   

        // int chr_len2=chr_tad_len[homology2/2];
        // int homology_start2=homology_tad_start[homology2];
        // int homology_end2=homology_start2+chr_len2;   
        // float dis;  
        // string t_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure_intra/test.txt";
        // ofstream test(t_file);        
        // for(int i=homology_start;i<homology_end;i++){
        //     for(int j=homology_start2;j<homology_end2;j++){
        //         // if(j>=homology_start && j<homology_end){continue;}
        //         test<<max(homology_dis_matrix[thread_num][min(i,j)][max(i,j)],float(1))<<'\t';
        //     }         
        //     test<<'\n';     
        // }
        // test.close();   
             

        // normalize_tads();
        // update_chrom_energies2();

              

        //输出homology_contact_matrix看一下结果
        string dis_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure_intra/dis.txt";
        ofstream distance(dis_file);
        for(int i=0;i<2958;i++){//2958
            for(int j=0;j<2958;j++){
                distance<<haploid_contact_matrix[i][j]<<'\t';
                // distance<<homology_dis_matrix[0][i][j]<<'\t';
            }
            distance<<'\n';
        }
        distance.close();                 


        //输出homology_dis_matrix看一下结果
        string c_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure_intra/contact.txt";
        ofstream contact(c_file);
        for(int i=0;i<500;i++){
            for(int j=0;j<500;j++){
                contact<<haploid_contact_matrix[i][j]<<'\t';
                // contact<<homology_contact_matrix[0][i][j]<<'\t';
            }
            contact<<'\n';
        }
        contact.close();              

        //输出chr_E看一下结果
        string e_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure_intra/e.txt";
        ofstream E(e_file);
        for(int i=0;i<500;i++){
            for(int j=0;j<500;j++){
                E<<chr_E[i][j]<<'\t';
            }
            E<<'\n';
        }
        E.close();   

        // get corr
        float cc=corr(haploid_contact_matrix);
        cout<<"intra_corr"<<'\n';
        cout<<"corr:"<<cc<<'\n';          

        // //输出homology_dis_matrix (染色体间contact)
        // string cont_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure/contact.txt";
        // ofstream cont(cont_file);

        // for(int i=0;i<22;i++){
        //     int chr_start=haploid_tad_start[i];
        //     int chr_end=chr_tad_len[i]+chr_start;
 
        //     for(int j=chr_start;j<chr_end;j++){
        //         for(int k=chr_start;k<chr_end;k++){
        //             haploid_contact_matrix[j][k]=0;
        //         }
        //     }
        // }

        // for(int i=0;i<2950;i++){//2958
        //     for(int j=0;j<2950;j++){
        //         cont<<haploid_contact_matrix[i][j]<<'\t';
        //         // distance<<homology_dis_matrix[0][i][j]<<'\t';
        //     }
        //     cont<<'\n';
        // }        
        // cont.close();                   

        //output polymer
        for(int i=0;i<5;i++){
            string p_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure_intra/polymer"+to_string(i)+".txt";
            ofstream P(p_file);
            for(int c=0;c<44;c++){
                int length=tads_center_chr[i][c].rows();
                MatrixXf &tmp = tads_center_chr[i][c];

                for (int l=0;l<length;l++){
                    P<<c<<'\t'<<tmp(l,0)<<'\t'<<tmp(l,1)<<'\t'<<tmp(l,2)<<'\n';
                }
            }
            P.close();              
        }

        //output polymer_bin
        for(int i=0;i<48;i++){
            string p_file="/data/home/txguo/code_final/reconstruct_genome_tad_structure_intra/polymer/polymer"+to_string(i)+".txt";
            ofstream P(p_file);
            for(int c=0;c<44;c++){
                int length=connected_structure_chr[i][c].rows();
                MatrixXi &tmp = connected_structure_chr[i][c];

                for (int l=0;l<length;l++){
                    P<<c<<'\t'<<tmp(l,0)<<'\t'<<tmp(l,1)<<'\t'<<tmp(l,2)<<'\n';
                }
            }
            P.close();              
        }        
                 

        // // get corr
        // float c=corr(haploid_contact_matrix);
        // cout<<"inter_corr"<<'\n';
        // cout<<"corr:"<<c<<'\n';  



    //     string time_interval=to_string(elapsed1.count());  
    //     corr_list.push_back({to_string(c),time_interval});

    //     //output 
    //     // string outfile=outfiledir+"/chr"+chr+"_rough.txt";
    //     ofstream fout(outfile_access);
    //     MatrixXf M = MatrixXf(tads_final_contacts_s)+reference_tads_contact_s.transpose();
    //     // MatrixXf M = tads_final_contacts_s+reference_tads_contact_s.transpose();
    //     for (int i=0;i<M.rows();++i){
    //         for(int j=0;j<M.rows();++j){
    //             fout << M(i,j) <<"\t";
    //         }
    //         fout<<'\n';
    //     }
    //     fout.close();     

        // //output polymer
        // out_polymer_filedir="/data/home/txguo/code_final/reconstruct_genome_tad_structure_intra/polymer";
        // for(auto l=0;l<number_of_threads;l++){
        //     get_final_structure(l,chr_homology,out_polymer_filedir);
        // }    
    }   


         
}