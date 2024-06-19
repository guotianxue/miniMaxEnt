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

int start_bin=0;
int end_bin=0;
int replicate_iter=-1;
//-------------外部参数-----------------

float learning_rate=0.4;//0.05
int radius=30;//the radius of cell nucleus
int number_of_threads = 48;//48
int update_steps =5;
long int burn_in_time = 10000;
long int mc_moves_start = 2000000;

int mask_bins=1;
vector<int> mask_bin_list;

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

void structure_initial(int start_bin,int end_bin,int model){//model=1 initial tads ; model=2 initial all chr
    //initial randomly
    int sitei,sitej;    

    cout<<start_bin<<"\t"<<end_bin<<"\n";

    //---------------------- initial ----------------------------

    tad_bin_num=end_bin-start_bin;
    pol_length = tad_bin_num*6;//monomers的数量

    //initial
    initial_oneTad_contact(start_bin,end_bin);

    //save refer
    string contacts_out=outfile+"/chr"+chr + "_refer/chr" + chr + "_" + to_string(start_bin)+ "_"+ to_string(end_bin-1)+ ".mtx" ;  
    SparseMatrix<float> contact=reference_one_tad_contact_s+SparseMatrix<float>(reference_one_tad_contact_s.transpose());
    Eigen::saveMarket(contact, contacts_out);       

    initial_mask(start_bin,end_bin);
    // reference_one_tad_contact_s=reference_one_tad_contact_s.cwiseProduct(mask_s);

    normalize_reference(reference_one_tad_contact_s) ; 
    initial_contact_mat(500);

    if (model==1){
        initial_tad_energy();
        for (int l = 0; l < number_of_threads; l++) {
            initial_polymer(polymer[l],position2site_index[l],site_index2position[l],site_available[l],pol_length,l); 
        }         
    }
    if (model==2){
        //get_tad data
        initial_tad_energy(); 
        cout<<"initial_e"<<'\n';
        for (int l = 0; l < number_of_threads; l++) {
            // string filename1="/data/home/txguo/data_use/maxEnt/whole_chr/rough_big_structure/2/"+to_string(l)+".txt";
            string filename1=inputfiledir+"/chr"+chr+"_"+to_string(l)+".txt";
            initial_polymer_byfile(polymer[l],position2site_index[l],site_index2position[l],site_available[l],pol_length,l,filename1);             
        }            
    }   

    cout << "Done with burn in " << endl;
    return ;
}

//------------ get corr ---------
float mean(SparseMatrix<float> &Matrix){
    return Matrix.sum()/(tad_bin_num*tad_bin_num);
}

float corr2(SparseMatrix<float> Matrix1,SparseMatrix<float> Matrix2){

    int length=Matrix1.rows();
    float mean1=mean(Matrix1);float mean2=mean(Matrix2);
    MatrixXf mean1_mat=MatrixXf::Constant(length,length, mean1);
    MatrixXf mean2_mat=MatrixXf::Constant(length,length, mean2);

    Matrix1=Matrix1-mean1_mat;
    Matrix2=Matrix2-mean2_mat;    
    
    float r=((Matrix1.cwiseProduct(Matrix2)).sum())/sqrt((Matrix1.cwiseProduct(Matrix1)).sum() * (Matrix2.cwiseProduct(Matrix2)).sum());
    return r;
}

void reconstruct(){;

    // show parameter
    if (initial_tads){cout<<"intial_tads:"<<'\n';}
    else{
        cout<<"intial_all_chr:"<<'\n';     
    }

    cout<<"#######################"<<'\n';
    cout<<"parameter:"<<'\n';
    cout<<"chr:"<<chr<<'\n';
    cout<<"radius:"<<'\t'<<radius<<'\n';
    cout<<"number_of_threads:"<<'\t'<<number_of_threads<<'\n';
    cout<<"update_steps:"<<'\t'<<update_steps<<'\n';
    // cout<<"burn_in_time:"<<'\t'<<burn_in_time<<'\n';
    cout<<"mc_moves_start:"<<'\t'<<mc_moves_start<<'\n';

    cout<<"neighbor_dis:"<<'\t'<<neighbor_dis<<'\n';

    initial_refer(chr);

    if(end_bin!=chr_index[chr]){
        end_bin=end_bin+1;//最后的一个bin也重建
    }        

    structure_initial(start_bin,end_bin,1);

    auto start = chrono::high_resolution_clock::now();

    for (int n =0; n<update_steps;n++) { //do iterative update scheme

        //------------------- initial ---------------------
        cout<<"---------------"<<'\n';
        // reduce_count=0,add_count=0,kink_count=0,loop_count=0,skip_count=0;

        cout<<"iter_step:"<< n<<"\n";
        mc_moves = mc_moves_start*sqrt(n+10)/sqrt(10); //number of MC moves grows with each iteration (implicitly converted to long int)

        tad_final_contacts_s=tad_final_contacts_s*0;
        //run forward simulation
        for (auto l = 0; l < number_of_threads; l++) {
            unordered_map<pair<int, int>, float, pair_hash>().swap(one_tad_contact[l]);
            unordered_map<int, vector<int>>().swap(tad_neighbor_contact[l]);
            unordered_map<pair<int, int>, NeighborDefaulted, pair_hash>().swap(neighbor_contact_data[l]);

            unordered_map<pair<int, int>, float, pair_hash>().swap(tad_total_contacts_neighbor[l]);
            threads[l] = thread(run, l, mc_moves);   

        }
        for (auto &&l : threads) {
            l.join();
        }

        // read in remaining contacts at the end of forward simulation
        for (auto l = 0; l < number_of_threads; l++) {
            for (auto elem : one_tad_contact[l]) { //add the contacts remaining at the end of the simulation
                // tad_total_contacts_s[l](elem.first.first/6,elem.first.second/6) += (mc_moves - elem.second);   
                tad_total_contacts_s[l].coeffRef(elem.first.first/6,elem.first.second/6) += (mc_moves - elem.second);
                // one_tad_contact[l][elem.first] = 0;
            }
        }

        // // read in remaining neighbor contacts at the end of forward simulation
        // for (auto l = 0; l < number_of_threads; l++) {
        //     for (auto elem : tad_neighbor_contact[l]) { //add the contacts remaining at the end of the simulation
        //         for (int k=0;k<elem.second.size();k++){
        //             int data=elem.second[k];
        //             tad_total_contacts_s[l].coeffRef(elem.first.first/6,elem.first.second/6) += float(mc_moves - neighbor_contact_m_mat[l](elem.first/6,data/6))*pow(neighbor_contact_dis_mat[l](elem.first/6,data/6),-3);                    
        //         }                        
        //     }
        // }

        // //get tad_total_contacts_neighbor
        // MatrixXf tad_final_contacts_neighbor=MatrixXf::Zero(tad_bin_num,tad_bin_num);
        // for(int l=0;l<number_of_threads;l++){
        //     for (auto map:tad_total_contacts_neighbor[l]){
        //         // map.first;
        //         tad_final_contacts_neighbor(map.first.first,map.first.second)=map.second;
        //     }
        // }
    

        //add up contacts from threads
        
        for (int l = 0; l <number_of_threads; l++) {
            tad_final_contacts_s += tad_total_contacts_s[l];  
        }

        // tad_final_contacts_s+=tad_final_contacts_neighbor.sparseView();
        tad_final_contacts_s=tad_final_contacts_s.cwiseProduct(mask_s);   

        auto time0 = chrono::high_resolution_clock::now();
        chrono::duration<float> elapsed0 = time0 - start;
        cout<<"run_time:"<<elapsed0.count() << " seconds\n";            

        //output contact row
        string contact="./e.txt";
        ofstream Cout(contact);  
        Cout<<MatrixXf(tad_final_contacts_s);
        Cout.close();                
        
        //normalize contact frequencies
        normalize();           

        //update energies
        update_energies(reference_one_tad_contact_s);                          

        auto time1 = chrono::high_resolution_clock::now();
        chrono::duration<float> elapsed1 = time1 - start;
        cout<<"run_time:"<<elapsed1.count() << " seconds\n";
        string time_interval=to_string(elapsed1.count());            

        float corr=corr2((tad_final_contacts_s+SparseMatrix<float>(tad_final_contacts_s.transpose())),(reference_one_tad_contact_s+SparseMatrix<float>(reference_one_tad_contact_s.transpose())));
        cout<<"corr:"<<corr<<'\n';     

        // string test="./test.txt";
        // ofstream Sout(test);

        // Sout<<MatrixXf(tad_final_contacts_s)+MatrixXf(reference_one_tad_contact_s).transpose()<<'\n';        
        // // Sout<<MatrixXf(tad_final_contacts_s);  
        // Sout.close() ;   


        // string polymer_out="./polymer.txt";
        // ofstream Pout(polymer_out);
        // int pol_index=0;
        // int l=0;
        // for (int i=0;i<site_available[l].size();i++){ 
        //     if (site_available[l][i]==1 ){

        //         Pout << polymer[l][pol_index][0] <<"\t"<< polymer[l][pol_index][1] <<"\t"<< polymer[l][pol_index][2] <<"\t"<< i <<"\n";
        //         pol_index+=1;
        //     }
        // }
        // Pout.close();     

        string contacts_out="/data/home/txguo/code_final/reconstruct_tads_pre/final.txt";
        ofstream Sout(contacts_out);
        // Sout<<(tad_final_contacts_s+SparseMatrix<float>(reference_one_tad_contact_s.transpose())).block(100,100,900,900);             
        Sout<<(tad_final_contacts_s+SparseMatrix<float>(reference_one_tad_contact_s.transpose())); 
        // Sout<<MatrixXf(reference_one_tad_contact_s);
        Sout.close();              
    }    

    if(end_bin!=chr_index[chr]){
        end_bin=end_bin-1;//方便命名
    }    

    //output contact
    string contacts_out=outfile+"/"+chr + "/chr" + chr + "_" + to_string(start_bin)+ "_"+ to_string(end_bin)+ ".mtx" ;  
    SparseMatrix<float> contact=tad_final_contacts_s+SparseMatrix<float>(reference_one_tad_contact_s.transpose());
    Eigen::saveMarket(contact, contacts_out);      

    //output polymer
    int pol_index=0;
    for(int l=0;l<number_of_threads;l++){
        string polymer_out=outfile+"/polymer/"+chr +"/chr"+chr+"_"+ to_string(start_bin)+ "_"+ to_string(end_bin)+ "_"+to_string(l)+".txt" ;
        ofstream Pout(polymer_out);     
        pol_index=0;   
        for (int i=0;i<site_available[l].size();i++){ 
            if (site_available[l][i]==1 ){

                Pout << polymer[l][pol_index][0] <<"\t"<< polymer[l][pol_index][1] <<"\t"<< polymer[l][pol_index][2] <<"\t"<< i <<"\n";
                pol_index+=1;
            }
        }
        Pout.close();  
    }
   

    // // //output polymer
    // string polymer_out="./polymer.txt";
    // // for(int l=0;l<number_of_threads;l++){
    // //     if (not initial_tads){
    // //         polymer_out = outfile+ "/polymer/chr" + chr +"_" + to_string(l) + ".txt" ;  
    // //     }
    // //     else{polymer_out = outfile+ "/polymer/chr" + chr + "_" + to_string(start_bin)+ "_"+ to_string(end_bin)  +"_" + to_string(l) + ".txt" ;  }

     
}


