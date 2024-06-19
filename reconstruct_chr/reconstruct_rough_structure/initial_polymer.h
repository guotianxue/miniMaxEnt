#ifndef Inverse_ploymer_Initialize_h
#define Inverse_ploymer_Initialize_h


#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <string>
#include "global.h"
#include "Energy_changes_liner.h"
using namespace std;


int tad_bin_num;
int pol_length;//monomers的数量

static mt19937_64 gen(time(0));//mt19937是伪随机数产生器，用于产生高性能的随机数

//=====================tads reconstruct initial=====================
Eigen::MatrixXf reference_tads_contact;//intial reference
Eigen::MatrixXd refer_contact;
Eigen::MatrixXf Interaction_E_tads;
unordered_map<int, int>  tad_bin_dict;
vector<vector<int>> tads_bound;
Eigen::MatrixXi tad_match;//(thread_nums,tads_bound.size());
Eigen::VectorXi threads_m(number_of_threads);

//initial contact and locations
vector<unordered_map<pair<int, int>, float, pair_hash>> tads_contact(number_of_threads);
vector<unordered_map<vector<int>, vector<int>, vec_hash>> tads_location(number_of_threads);
vector<vector<Eigen::MatrixXi>> tads(number_of_threads);
vector<vector<Eigen::MatrixXi>> connected_tads(number_of_threads);
vector<Eigen::MatrixXi> connected_structure(number_of_threads);
vector<Eigen::MatrixXf> tads_total_contacts(number_of_threads);
Eigen::MatrixXf tads_final_contacts;

Eigen::MatrixXd final_contacts;

//tad_center_initial
vector<Eigen::MatrixXf> tads_center(number_of_threads);
vector<Eigen::MatrixXf> tads_center_m(number_of_threads);
vector<Eigen::MatrixXf> tads_center_dis(number_of_threads);

//mutex use
vector<int> is_move_active(number_of_threads);//is move working
int tad_change_flag=0;
int tad_change_thread=-1;



//------------------------------polymer初始化------------------------------------------------
bool check_nucleus(int x,int y,int z){
    return (pow(x,2)+pow(y,2)+pow(z,2)<=pow(radius,2));
}

// //============================= tads calculate ============================

void initial_rough_refer(string chr,int tad_num){
    //get a rough structure (a point is a tad)    
    Eigen::MatrixXf contact=MatrixXf::Zero(tad_num,tad_num);
    
    string refer_tad_file="/data/home/txguo/data_use/maxEnt/whole_chr/tad_contact/chr"+(chr)+"_tad_contact.txt";

    ifstream tad_refer;
    tad_refer.open(refer_tad_file);

    int sitei;
    int sitej;
    double contact_num;

    while(tad_refer >> sitei >> sitej >> contact_num){
        if(sitei!=sitej){
            contact(sitei,sitej)=contact_num;
        }
    }   
    
    tad_refer.close();     
    reference_tads_contact=contact.block(0,0,tads_bound.size(),tads_bound.size());

}

// //------------------------deal with tad_level data-------------------------

void connect_tad(Eigen::MatrixXi first_point,Eigen::MatrixXi &tad,Eigen::MatrixXi &connected_tad){
    vector<Eigen::Matrix<int,1,3>> rotate={{1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};    
    Eigen::MatrixXi margin=first_point.replicate(tad.rows(), 1);
    Eigen::MatrixXi rotate_mat;
    //rotate
    int dir=unipath(gen);
    rotate_mat=rotate[dir].replicate(tad.rows(), 1);
    connected_tad=tad.cwiseProduct(rotate_mat) +margin;   //dot                    
}
void initial_tads(int thread_num){
    //initial tad center
    tads_center[thread_num]=Eigen::MatrixXf::Zero(tads_bound.size(),3);
    Eigen::MatrixXi  x;
    Eigen::MatrixXi  y;
    Eigen::MatrixXi  z;
    int first=tads_bound[0][0];
    int start;int end;float x_c;float y_c;float z_c;
    for (int i=0;i<tads_bound.size();i++){
        start=tads_bound[i][0]-first;
        end=tads_bound[i][1]-first;
        x=connected_structure[thread_num].middleRows(start,end-start).col(0);
        y=connected_structure[thread_num].middleRows(start,end-start).col(1);
        z=connected_structure[thread_num].middleRows(start,end-start).col(2);
        x_c=float(x.sum())/(end-start);
        y_c=float(y.sum())/(end-start);
        z_c=float(z.sum())/(end-start);
        tads_center[thread_num].row(i)<<x_c,y_c,z_c;
        }   

    //initial the dis to the center
    tads_center_dis[thread_num]= Eigen::MatrixXf::Zero(tads_bound.size(),tads_bound.size()); 
    Eigen::MatrixXf margin;
    for (int i=0;i<tads_bound.size();++i){
        for(int j=i+1;j<tads_bound.size();++j){
            margin=tads_center[thread_num].row(j)-tads_center[thread_num].row(i);
            tads_center_dis[thread_num](i,j)=sqrt((margin.cwiseProduct(margin)).sum());    
            tads_center_dis[thread_num](j,i)=tads_center_dis[thread_num](i,j);            
        }
    }

    //initial tads total contacts
    int length=tads_bound[tads_bound.size()-1][1]+1-tads_bound[0][0];
    tads_total_contacts[thread_num]=Eigen::MatrixXf::Zero(tads_bound.size(),tads_bound.size()); 

}

void initial_tads_data(){
    //initial match
    Eigen::MatrixXi thread_match(number_of_threads,1);
    int start,end;

    //initial tad_match 
    for(int i=0;i<number_of_threads;i++){
        thread_match(i,0)=i;
    }       
    tad_match=thread_match.replicate(1, tads_bound.size()); 

    //initial m of every threads
    threads_m=Eigen::VectorXi::Zero(number_of_threads);

    //initial tad dict
    for (int i=0;i<tads_bound.size();i++){
        start=tads_bound[i][0];
        end=tads_bound[i][1];
        for (int j=start;j<end;j++){
            tad_bin_dict[j]=i;
        }    

    tad_bin_dict[tads_bound[tads_bound.size()-1][1]]=tads_bound.size()-1;


    }

    //initial tad energy
    int length=tads_bound[tads_bound.size()-1][1]+1-tads_bound[0][0];
    Interaction_E_tads=Eigen::MatrixXf::Zero(tads_bound.size(),tads_bound.size());
}

// //read tad bin from file
void initial_tads_byfile(int thread_num,vector<Eigen::MatrixXi> &tads,vector<Eigen::MatrixXi> &connected_tads,Eigen::MatrixXi &connected_structure,string files,vector<vector<int>> tads_bound){
    int x,y,z,loc;
    int x_ref,y_ref,z_ref;
    int start,end;
    //initial first point of first tad
    x=uniradius(gen);
    y=uniradius(gen);
    z=uniradius(gen);
    while(check_nucleus(x,y,z)==0){
        x=uniradius(gen);
        y=uniradius(gen);
        z=uniradius(gen);       
    }    
    Eigen::MatrixXi tad;
    Eigen::MatrixXi connected_tad;
    Eigen::MatrixXi first_point(1,3);
    first_point<<x,y,z;
    // save tad and get connected tad
    ifstream initial_polymer;
    for(auto bound: tads_bound ){
        int i=0;
        start=bound[0];
        end=bound[1];
        // cout<<start<<'\t'<<end<<'\n';
        string tad_file=files+"/polymer_"+to_string(start)+"_"+to_string(end)+"_"+to_string(thread_num)+".txt";
        Eigen::MatrixXi tad=Eigen::MatrixXi::Zero(end-start+1,3);
        Eigen::MatrixXi connected_tad=Eigen::MatrixXi::Zero(end-start+1,3);

        initial_polymer.open(tad_file);
        while(initial_polymer>> x >> y >> z >> loc){
            if (i==0){
                x_ref=x;
                y_ref=y;
                z_ref=z;
                }
            if (loc%6==0){
                //set the first point is zero
                tad(i,0)=x-x_ref;
                tad(i,1)=y-y_ref;
                tad(i,2)=z-z_ref;
                i+=1;

            }
        }
        initial_polymer.close();
        //-----connect tads and initial contacts-----
        connect_tad(first_point,tad,connected_tad);
        tads.push_back(tad);
        connected_tads.push_back(connected_tad);
        first_point=connected_tad(connected_tad.rows()-1,Eigen::all);
    }

    //initial big structure 
    vector<int> monomer;
    int first= tads_bound[0][0];
    // cout<<end<<'\n';
    connected_structure.resize(end - first +1,3);

    for(int i=0;i<tads_bound.size();i++){
        start=tads_bound[i][0];
        end=tads_bound[i][1];
        connected_structure.block(start-first,0,end-start+1,3)=connected_tads[i];  
        Eigen::MatrixXi &tad_use=connected_tads[i];

        //update location 和 contact 
        for (int j=start;j<end;j++){
            monomer={tad_use(j-start,0),tad_use(j-start,1),tad_use(j-start,2)};
            if (tads_location[thread_num].find(monomer)!=tads_location[thread_num].end()){
                for (auto elem : tads_location[thread_num][monomer]){
                    if (elem<start)
                        tads_contact[thread_num][{min(elem,j),max(elem,j)}]=0;//这个0是用来记录后面第几次move  
                }
                tads_location[thread_num][monomer].push_back(j);//三维位置重合的点加入Location
            }
            else { tads_location[thread_num][monomer] = {j}; }//如果monomer不在location中，即第一次出现，加入该monomer            
        }
    }

    //initial tad data
    initial_tads(thread_num);

}

//===========================================================
void get_tadset(string filepath, vector<vector<int>> &tads_bound){
    ifstream path;
    path.open(filepath);
    int sitei,sitej;
    int index=0;
       
    while(path >> sitei >> sitej ){     
        tads_bound.push_back({sitei,sitej});
        index+=1;
    }          
    path.close();       
}


#endif