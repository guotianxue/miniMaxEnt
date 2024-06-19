#ifndef Inverse_ploymer_Initialize_h
#define Inverse_ploymer_Initialize_h


#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <string>
#include "global.h"
#include "Energy_changes_liner.h"
using namespace std;
using namespace Eigen;


int chr_num=22;
int chr_homology_num=44;
int tads_num=0;

vector<int> connect_dir={0,1,2,3,4,5,6,7};
vector<vector<MatrixXi>>  dir_match(number_of_threads,vector<MatrixXi>(chr_homology_num));
//对应rotate={{1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}}

unordered_map<int, int>  chr_len;//chr_index:chr_length
unordered_map<string, int>  chr_index_dict;//chr_name:chr_index
unordered_map<int, string> chr_name_dict;

vector<int> chr_tad_len(chr_num,0);
vector<int> haploid_tad_start(chr_num);
vector<int> homology_tad_start(chr_homology_num);

vector<unordered_map<int, vector<int>>>  chr_tad2bound(chr_num); //chr_tad_index:{tad_bin_start,tad_bin_end}
unordered_map<int, vector<int>>  genome_tad2bound;//genome_tad_index:{tad_bin_start,tad_bin_end}这个不含同源染色体
unordered_map<int, vector<int>>  homology_tad2chr_tad; //homology_tad_index:{chr_homology,chr_tad_index}
// unordered_map<vector<int>, int, vec_hash> chr_tad2homology_tad; //{chr_homology,chr_tad_index}:homology_tad_index
vector<vector<int>> chr_tad2homology_tad(chr_homology_num);

vector<vector<vector<float>>> homology_dis_matrix(number_of_threads);
vector<vector<vector<float>>> homology_dis_matrix_change(number_of_threads);

vector<vector<int>> homology_m(number_of_threads);//用完要更新成全0
vector<int> thread_m(number_of_threads,0);//用完要更新成全0

vector<vector<vector<float>>> homology_contact_matrix(number_of_threads);
vector<vector<float>> haploid_contact_matrix;

vector<vector<int>> chr_contact_matrix;

vector<vector<float>> chr_E;

static mt19937_64 gen(time(0));//mt19937是伪随机数产生器，用于产生高性能的随机数

//=====================tads reconstruct initial=====================
// Eigen::MatrixXf reference_tads_contact;//intial reference
Eigen::MatrixXd refer_contact;
Eigen::MatrixXf Interaction_E_tads;
Eigen::MatrixXd Interaction_E;
unordered_map<int, int>  tad_bin_dict;
vector<vector<int>> tads_bound;
vector<vector<int>> tads_bound_genome;
vector<vector<int>> tads_bound_chr;
vector<Eigen::MatrixXi> tad_match_chr(chr_homology_num);//[chr_homology](thread_nums,tads_bound.size());

//initial sparse matrix
SparseMatrix<float> reference_genome_tads_contact_s;
SparseMatrix<float> reference_tads_contact_s;
vector<vector<float>> reference_genome_tads_contact;


//initial contact and locations
vector<unordered_map<pair<int, int>, float, pair_hash>> tads_contact(number_of_threads);
vector<unordered_map<vector<int>, vector<int>, vec_hash>> tads_location(number_of_threads);

vector<vector<vector<Eigen::MatrixXi>>> tads_chr(number_of_threads,vector<vector<Eigen::MatrixXi>>(chr_homology_num));
vector<vector<vector<Eigen::MatrixXi>>> connected_tads_chr(number_of_threads);//connected_tads[thread_num][chri][index]
vector<vector<Eigen::MatrixXi>> connected_structure_chr(number_of_threads,vector<Eigen::MatrixXi>(chr_homology_num));
vector<Eigen::MatrixXi> connected_structure_genome(number_of_threads);

vector<vector<Eigen::MatrixXi>> connected_tads(number_of_threads);
vector<Eigen::MatrixXi> connected_structure(number_of_threads);

vector<SparseMatrix<float>> tads_total_contacts_s(number_of_threads);

SparseMatrix<float> tads_final_contacts_s;

Eigen::MatrixXf mask;
SparseMatrix<float> mask_s;

Eigen::MatrixXd final_contacts;

//chrom_point_initial
vector<vector<MatrixXi>> chrom_point(number_of_threads);

//tad_center_initial
vector<Eigen::MatrixXf> tads_center(number_of_threads);
vector<Eigen::MatrixXf> tads_center_m(number_of_threads);
vector<Eigen::MatrixXf> tads_center_dis(number_of_threads);
vector<vector<Eigen::MatrixXf>> tads_center_chr(number_of_threads,vector<Eigen::MatrixXf>(chr_homology_num));

//mutex use
vector<int> is_move_active(number_of_threads);//is move working
int tad_change_flag=0;
int tad_change_thread=-1;



//------------------------------polymer初始化------------------------------------------------
bool check_nucleus(int x,int y,int z){
    return (pow(x,2)+pow(y,2)+pow(z,2)<=pow(radius,2));
}



//===================data_initial===============================

void initial_chrom_point(){
    string chrom_initial_dir="/data/home/txguo/data_use/maxEnt/whole_chr/chrom_initial/";
    string chrom_initial_file;
    int x,y,z,index=0;
    MatrixXi use(1,3);

    for(int i=0;i<number_of_threads;++i){

        chrom_initial_file=chrom_initial_dir+"chrom_initial_"+to_string(i)+".txt";
        ifstream chrom_initial;
        chrom_initial.open(chrom_initial_file);

        while(chrom_initial >> x >> y >> z){
            use<<x,y,z;
            chrom_point[i].push_back(use);
        } 
    }

}

void initial_chrom_point_zero(){
    int x,y,z,index=0;
    MatrixXi use(1,3);
    use<<0,0,0;

    for(int i=0;i<number_of_threads;++i){
        for(int j=0;j<chr_homology_num;j++){
            chrom_point[i].push_back(use);
        }
    }
}

void initial_tad_contact_refer(){
    //get a rough structure (a point is a tad)    

    string genome_tad_contact_file="/data/home/txguo/data_use/maxEnt/whole_chr/tad_contact/genome_tad_contact.mtx";
  
    Eigen::loadMarket(reference_genome_tads_contact_s, genome_tad_contact_file);      
}

void cut_tads_refer(int chr_homology){
    
    int tad_start=haploid_tad_start[chr_homology/2];;
    int chr_len=chr_tad_len[chr_homology/2];
    reference_tads_contact_s=reference_genome_tads_contact_s.block(tad_start,tad_start,chr_len,chr_len);
    cout<<reference_tads_contact_s.cols()<<'\t'<<chr_len<<'\n';

}

void get_tadset(vector<vector<int>> &tads_bound){

    string filepath="/data/home/txguo/data_use/maxEnt/whole_chr/divide_tad/tadlevel_divide.txt";  
    ifstream path;
    path.open(filepath);
    string chr="chr",chr_pre="";
    int sitei,sitej,sitem,siten,chr_index=-1,start_index=0;
    int index=0;
    
    //initial chr_index_dict
    while(path >> sitei >> sitej >> chr >> sitem >> siten){   
        if(chr_pre!=chr) {
            chr_index++;
            chr_index_dict[chr]=chr_index;
            chr_name_dict[chr_index]=chr;
            haploid_tad_start[chr_index]=start_index;

        }

        chr_tad2bound[chr_index][start_index-haploid_tad_start[chr_index]]={sitem,siten};
        genome_tad2bound[start_index]={sitei,sitej};      
        start_index++;


        chr_pre=chr;
        chr_tad_len[chr_index]+=1;
        tads_bound_genome.push_back({sitei,sitej});
        tads_bound_chr.push_back({chr_index,sitem,siten});
        
    }          
    path.close();  

    //initial chr_len
    int chr_length;
    string chr_size_file="/data/home/txguo/data_use/maxEnt/reference/reference_use/hg19.chrom_20k.sizes";
    
    ifstream chr_size;
    chr_size.open(chr_size_file);
    int i=0;
    while(chr_size >> chr>> chr_length){
        chr_len[chr_index_dict[chr]]=chr_length;
        i++;
        if(i==chr_num){
            break;
        }
    }
    chr_size.close();  

    //initial
    int homology_tad_index=0;
    for(int chr_homology=0;chr_homology<chr_homology_num;chr_homology++){
        chr_index=chr_homology/2;
        int chr_tad_num=chr_tad_len[chr_index];
        
        for(int tad_index=0;tad_index<chr_tad_num;tad_index++){

            if(tad_index==0){
                homology_tad_start[chr_homology]=homology_tad_index;
            }

            homology_tad2chr_tad[homology_tad_index]={chr_homology,tad_index};
            chr_tad2homology_tad[chr_homology].push_back(homology_tad_index) ;

            homology_tad_index++;
        }
    }

    //initial tads_num
    tads_num=homology_tad2chr_tad.size();


}

// //------------------------deal with tad_level data-------------------------
void move_tad(Eigen::MatrixXi move_point,Eigen::MatrixXi &tad,int multiple){

    Eigen::MatrixXi margin=move_point.replicate(tad.rows(), 1)*multiple;
    tad=tad + margin;                    
}

void move_center(Eigen::MatrixXf move_point,Eigen::MatrixXf &tad,float multiple){

    Eigen::MatrixXf margin=move_point.replicate(tad.rows(), 1)*multiple;
    tad=tad + margin;                    
}


void move_connect_tad(Eigen::MatrixXi move_point,Eigen::MatrixXi &tad,Eigen::MatrixXi &connected_tad,int dir){
    vector<Eigen::Matrix<int,1,3>> rotate={{1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};    
    Eigen::MatrixXi margin=move_point.replicate(tad.rows(), 1);
    Eigen::MatrixXi rotate_mat;
    //rotate
    rotate_mat=rotate[dir].replicate(tad.rows(), 1);
    connected_tad=tad.cwiseProduct(rotate_mat) +margin;   //dot                  
}

float get_connect_dis(Eigen::MatrixXi move_point,Eigen::MatrixXi &tad,int dir){

    MatrixXi x,y,z,connected_tad;
    float x_c,y_c,z_c,dis;
    int length=tad.rows();

    vector<Eigen::Matrix<int,1,3>> rotate={{1,1,1},{1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,1},{-1,-1,1},{-1,1,-1},{-1,-1,-1}};    
    Eigen::MatrixXi rotate_mat;
    //rotate
    rotate_mat=rotate[dir].replicate(tad.rows(), 1);    
    connected_tad=tad.cwiseProduct(rotate_mat);
    x=connected_tad.col(0);
    y=connected_tad.col(1);
    z=connected_tad.col(2);
    x_c=float(x.sum())/(length);
    y_c=float(y.sum())/(length);
    z_c=float(z.sum())/(length);
    dis=sqrt((x_c+move_point(0,0))*(x_c+move_point(0,0))+(y_c+move_point(0,1))*(y_c+move_point(0,1))+(z_c+move_point(0,2))*(z_c+move_point(0,2)));

    return dis;
}


void initial_tads_center(int thread_num){
    //initial tad center
    int homology_chr_index;
    Eigen::MatrixXi  x,y,z;
    int start,end;
    float x_c;float y_c;float z_c;    

    for(int chr_index=0;chr_index<chr_num;chr_index++){
        for(int homology=0;homology<2;homology++){
            homology_chr_index=chr_index*2+homology;
            tads_center_chr[thread_num][homology_chr_index]=Eigen::MatrixXf::Zero(chr_tad_len[chr_index],3);           
        }
    
    }

    int chr_index;
    for (int i=0;i<tads_bound_chr.size();i++){
        for(int homology=0;homology<2;homology++){
            chr_index=tads_bound_chr[i][0];
            start=tads_bound_chr[i][1];
            end=tads_bound_chr[i][2];
            homology_chr_index=chr_index*2+homology;

            x=connected_structure_chr[thread_num][homology_chr_index].middleRows(start,end-start).col(0);
            y=connected_structure_chr[thread_num][homology_chr_index].middleRows(start,end-start).col(1);
            z=connected_structure_chr[thread_num][homology_chr_index].middleRows(start,end-start).col(2);
            x_c=float(x.sum())/(end-start);
            y_c=float(y.sum())/(end-start);
            z_c=float(z.sum())/(end-start);

            tads_center_chr[thread_num][homology_chr_index].row(i-haploid_tad_start[chr_index])<<x_c,y_c,z_c;  

        }
    }    

    // move center point to  chrom point  
    MatrixXi move_coordinate,chrom_coordinate;
    int mid_tad_index,mid_point_index;     
    for(int chr_index=0;chr_index<chr_num;chr_index++){
        // int a=chr_tad_len[chr_index];    
        // cout<<chr_index<<'\t'<<a<<'\t'<<length<<'\n';
        for(int homology=0;homology<2;homology++){
            homology_chr_index=chr_index*2+homology;   
            int length=tads_center_chr[thread_num][homology_chr_index].rows();
            x_c=tads_center_chr[thread_num][homology_chr_index].col(0).sum()/(length);
            y_c=tads_center_chr[thread_num][homology_chr_index].col(1).sum()/(length);
            z_c=tads_center_chr[thread_num][homology_chr_index].col(2).sum()/(length);

            // cout<<"before"<<homology_chr_index<<'\t'<<x_c<<'\t'<<y_c<<'\t'<<z_c<<'\n';

            MatrixXi mid_point(1,3);
            mid_point<<x_c,y_c,z_c;

            //get chrom_point coordinate
            chrom_coordinate=chrom_point[thread_num][homology_chr_index];
            move_coordinate=chrom_coordinate-mid_point;   
            move_tad(move_coordinate,connected_structure_chr[thread_num][homology_chr_index],1); 
            Matrix<float,1,3> move_center_coordinate(float(move_coordinate(0,0)),float(move_coordinate(0,1)),float(move_coordinate(0,2)));
            // cout<<move_center_coordinate<<'\n';
            // cout<<tads_center_chr[thread_num][homology_chr_index](0,0)<<'\t'<<tads_center_chr[thread_num][homology_chr_index](0,1)<<'\t'<<tads_center_chr[thread_num][homology_chr_index](0,2)<<'\n';            
            move_center(move_center_coordinate,tads_center_chr[thread_num][homology_chr_index],1); 
            // cout<<tads_center_chr[thread_num][homology_chr_index](0,0)<<'\t'<<tads_center_chr[thread_num][homology_chr_index](0,1)<<'\t'<<tads_center_chr[thread_num][homology_chr_index](0,2)<<'\n';
        }     

    }
}


void initial_tads_dis(int thread_num,int chr_homology){
    
    //initial dis_matrix and dis_matrix change
    int tads_num=chr_tad_len[chr_homology/2];
    vector<vector<float>> dis_matrix(tads_num,vector<float> (tads_num,0));
    homology_dis_matrix[thread_num]=dis_matrix;      

    // //initial the dis to the center
    int homology1,homology2,homology_index1,homology_index2;
    float dis;
    MatrixXf margin;

    for(int i=0;i<tads_num;i++){
        for(int j=i+1;j<tads_num;j++){
            margin=tads_center_chr[thread_num][chr_homology](i,all)-tads_center_chr[thread_num][chr_homology](j,all);
            homology_dis_matrix[thread_num][i][j] = sqrt(margin(0,0)*margin(0,0)+margin(0,1)*margin(0,1)+margin(0,2)*margin(0,2));
            
        }
    }       
    
    // for(int i=0;i<tads_num;i++){
    //     homology1=homology_tad2chr_tad[i][0];
    //     homology_index1=homology_tad2chr_tad[i][1];

    //     for(int j=i+1;j<tads_num;j++){
    //         homology2=homology_tad2chr_tad[j][0];
    //         homology_index2=homology_tad2chr_tad[j][1];

    //         margin=tads_center_chr[thread_num][homology1](homology_index1,all)-tads_center_chr[thread_num][homology2](homology_index2,all);

    //         homology_dis_matrix[thread_num][i][j] = sqrt(margin(0,0)*margin(0,0)+margin(0,1)*margin(0,1)+margin(0,2)*margin(0,2));
            
    //     }
    // }          

}

//read tad bin from tad file
void initial_tads_byfile(int thread_num,vector<vector<Eigen::MatrixXi>> &tads_chr,vector<Eigen::MatrixXi> &connected_structure_chr,string files){
    int x,y,z,loc,i,dir;
    int start,end,start_index=0,chr_tad_index;
    int chr_index=-1,tad_index=-1,chr_length,chr_homology;
    string chr;


    Eigen::MatrixXi first_point(1,3);
    Eigen::MatrixXi tad;
    vector<Eigen::MatrixXi> connected_tad(2);
    first_point<<0,0,0;        
 

    // save tad and get connected tad  
    ifstream initial_polymer;

    for(int k=0;k<chr_homology_num;++k){
        chr_index=k/2;
        chr_length=chr_len[chr_index];
        connected_structure_chr[k].resize(chr_length,3);
        connected_structure_chr[k]=Eigen::MatrixXi::Zero(chr_length,3); 

        dir_match[thread_num][k]=Eigen::MatrixXi::Zero(1,chr_tad_len[chr_index]); 
    }

    for(auto bound: tads_bound_chr ){
        tad_index+=1;

        chr_index=bound[0];
        chr=chr_name_dict[chr_index];
  
        start=bound[1];
        end=bound[2];

        string tad_file=files+"/"+chr+"/chr"+chr+"_"+to_string(start)+"_"+to_string(end)+"_"+to_string(thread_num)+".txt";
        //initial
        if(haploid_tad_start[chr_index]+chr_tad_len[chr_index]==tad_index+1){
            // cout<<haploid_tad_start[chr_index]<<'\t'<<chr_tad_len[chr_index]<<'\t'<<tad_index<<'\n';
            end--;
        }
        tad=MatrixXi::Zero(end-start+1,3);

        i=0;                
        initial_polymer.open(tad_file);
        while(initial_polymer>> x >> y >> z >> loc){
            if (loc%6==0){
                //initial the first point of tad is {0,0,0}
                tad(i,0)=x;
                tad(i,1)=y;
                tad(i,2)=z;
                i+=1;
            }
        }

        initial_polymer.close();

        for(int homology=0;homology<2;homology++){//homology initial

            chr_homology=chr_index*2+homology;

            
            move_tad(tad(0,Eigen::all),tad,-1);

            tads_chr[chr_homology].push_back(tad);
            
            //-----connect tads and initial contacts-----

            if(haploid_tad_start[chr_index]==tad_index ){
                first_point<<0,0,0;                
                if(homology==0){
                    chr_tad_index=0;
                }        
                else   chr_tad_index++;      
            }
            else{
                chr_tad_index++;
                first_point=connected_tad[homology](connected_tad[homology].rows()-1,Eigen::all);
            }

            uniform_int_distribution<int> uniconnect(0,7);
            int dir=uniconnect(gen);
            dir_match[thread_num][chr_homology](0,chr_tad_index/2)=dir; 
            //确保结构小于半径
            float dis=get_connect_dis(first_point,tad,dir);
            float dis_change;
            if(dis>radius){
                for(int i=0;i<=7;i++){
                    dis_change=get_connect_dis(first_point,tad,i);
                    if(dis_change<dis){
                        dis=dis_change;
                        dir=i;
                    }
                }
            }            

            move_connect_tad(first_point,tad,connected_tad[homology],dir);   
            connected_structure_chr[chr_homology].block(start,0,end-start+1,3)=connected_tad[homology];              
 
        }
    }    

    //initial tad center and move the center point 
    initial_tads_center(thread_num);

    // //initial dis_matrix and dis_matrix change
    // vector<vector<float>> dis_matrix(tads_num,vector<float> (tads_num,0));
    // homology_dis_matrix[thread_num]=dis_matrix;      

}


void initial_tads_data(){

    //initial tad_match 
    MatrixXi thread_match(number_of_threads,1);
    for(int k=0;k<chr_homology_num;k++){
        int chr_index=k/2;
        for(int i=0;i<number_of_threads;i++){
            thread_match(i,0)=i;
        } 
        tad_match_chr[k]=thread_match.replicate(1, chr_tad_len[chr_index]);       
    }

    //initial tad energy
    chr_E=vector<vector<float>>(tads_bound_genome.size(),vector<float>(tads_bound_genome.size(),0));

    
}

#endif