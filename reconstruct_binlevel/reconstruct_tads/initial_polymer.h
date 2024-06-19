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
// #include "Energy_changes_liner.h"
using namespace std;
using namespace Eigen;


int tad_bin_num;
int pol_length;//monomers的数量
int chr_len;

static mt19937_64 gen(time(0));//mt19937是伪随机数产生器，用于产生高性能的随机数

vector<vector<vector<int>>> polymer(number_of_threads);
vector<vector<vector<int>>> polymer_bin(number_of_threads);

//record polymer site index (vector) ------> (site index) : (polymer position site)
vector<vector<int>> site_index2position(number_of_threads); //(number_of_threads,vector<int>(pol_length,0));

//record polymer site index (vector) ------> (polymer position site) : (site index)
vector<vector<int>> position2site_index(number_of_threads); //(number_of_threads,vector<int>(pol_length,-1));

//record site available flag     0 means no site   1 means have site
vector<vector<int>> site_available(number_of_threads); //(number_of_threads,vector<int>(pol_length,0));

//energy initial 
SparseMatrix<float> Interaction_E_tad_s;
Eigen::MatrixXf Interaction_E_tad;
unordered_map<pair<int, int>, float, pair_hash> Interaction_E_tad_1;
unordered_map<pair<int, int>, float, pair_hash>* volatile E_p=&Interaction_E_tad_1;

//initial reference
SparseMatrix<float> reference_contact_s; 
SparseMatrix<float> reference_one_tad_contact_s;
Eigen::MatrixXf refer_msk_contact;  
Eigen::MatrixXf  ref_final;

unordered_map<string,int> chr_index;

//initial mask
SparseMatrix<float> mask_s;

//initial  multiple matrix
Eigen::MatrixXf matrix_multiple;
Eigen::MatrixXf deal_multiple;

//initial contact and locations
vector<unordered_map<pair<int, int>, float, pair_hash>> one_tad_contact(number_of_threads);
vector<unordered_map<vector<int>, vector<int>, vec_hash>> one_tad_locations(number_of_threads);

vector<SparseMatrix<float>> tad_total_contacts_s(number_of_threads);
vector<SparseMatrix<float>> tad_neighbor_contacts_s(number_of_threads);
SparseMatrix<float> tad_final_contacts_s;

//initial neighbor
vector<unordered_map<pair<int, int>, float, pair_hash>> tad_total_contacts_neighbor(number_of_threads);
vector<unordered_map<pair<int, int>, float, pair_hash>> neighbor_contact_dis_hash(number_of_threads);
vector<unordered_map<pair<int, int>, int, pair_hash>> neighbor_contact_m_hash(number_of_threads);
vector<unordered_map<pair<int, int>, NeighborDefaulted, pair_hash>> neighbor_contact_data(number_of_threads);

vector<unordered_map<int, vector<int>>> tad_neighbor_contact(number_of_threads);
vector<vector<int>> tad_neighbor_contact_new(number_of_threads);

vector<unordered_map<pair<int, int>, float, pair_hash>> neighbor_contact_dis(number_of_threads);

unordered_map<pair<int, int>, float, pair_hash> refer_hash;

//------------------------------polymer初始化------------------------------------------------
bool check_nucleus(int x,int y,int z){
    return (pow(x,2)+pow(y,2)+pow(z,2)<=pow(radius,2));
}

//Normalizes model contact frequencies to allow for comparison with experimental data
void normalize() {

    float sum = tad_final_contacts_s.sum();
    float check1=0;
    float multiply= float(tad_bin_num)/(2*sum);

    tad_final_contacts_s*= multiply;
    cout<<"nor_check"<<"\t"<<tad_final_contacts_s.sum()<<"\n";

    // float sum_refer=0;
    // float sum=0;
    // for (int i=0;i<tad_bin_num-2;i++){
    //     sum_refer+=reference_one_tad_contact_s.coeffRef(i,i+2);
    //     sum+=tad_final_contacts_s.coeffRef(i,i+2);
    // }
    // cout<<sum_refer<<'\t'<<sum<<'\n';
    // tad_final_contacts_s*= float(sum_refer)/(sum);
    
}

void normalize_reference(SparseMatrix<float> &refer) {

    float sum=refer.sum();
    refer*=tad_bin_num/(2*sum);
    
    cout<<"refer_check"<<"\t"<<refer.sum()<<"\n";
}

void initial_polymer(vector<vector<int>> &polymer,vector<int> &position2site_index,vector<int> &site_index2position,vector<int> &site_available,int pol_length,int thread_num){//随机游走初始化

    site_available.resize(pol_length);
    position2site_index.resize(pol_length);
    site_index2position.resize(pol_length);

    fill(site_available.begin(), site_available.end(), 0);
    fill(position2site_index.begin(), position2site_index.end(), 0);
    fill(site_index2position.begin(), site_index2position.end(), 0);
    //site_available initial and site_index2position initial 
    for (int i=0;i<tad_bin_num;++i){
        int tmp=6*i;
        int tmp1=4*i;
        //site available 100111循环
        site_available[tmp]=1;
        site_available[tmp+3]=1;
        site_available[tmp+4]=1;
        site_available[tmp+5]=1; 

        position2site_index[tmp1]=tmp;
        position2site_index[tmp1+1]=tmp+3;
        position2site_index[tmp1+2]=tmp+4;  
        position2site_index[tmp1+3]=tmp+5;   
                      

        if (i==0) {
            // site_index2position[0]=0;
            for (int j=1;j<6;j++){
                site_index2position[tmp+j]=site_index2position[tmp+j-1]+site_available[tmp+j];
            }            
        }
        else{
            for (int j=0;j<6;j++){
                site_index2position[tmp+j]=site_index2position[tmp+j-1]+site_available[tmp+j];
            }
        }
    }

    vector<vector<int>>().swap(polymer);
    vector<vector<int>>().swap(polymer_bin[thread_num]);
    unordered_map<vector<int>, vector<int>, vec_hash>().swap(one_tad_locations[thread_num]);
    unordered_map<pair<int, int>, float, pair_hash>().swap(one_tad_contact[thread_num]);


    //initial first point
    int x;int y;int z;
    x=uniradius(gen);
    y=uniradius(gen);
    z=uniradius(gen);
    while(check_nucleus(x,y,z)==0){
        x=uniradius(gen);
        y=uniradius(gen);
        z=uniradius(gen);       
    }  
    // x=0;y=0;z=0;
    vector<int> monomer={x,y,z}; //第一个点从原点出发
    polymer.push_back(monomer);
    polymer_bin[thread_num].push_back({x,y,z});
    one_tad_locations[thread_num][monomer] = {0};
        
    int k=0; int x_tmp;int y_tmp;int z_tmp;int path;
    for (int i=1;i<tad_bin_num*4;i++){//initial state: each gene site have 3 site interval
        while(1){
            x_tmp=x;y_tmp=y;z_tmp=z;
            path=unipath(gen);  
            switch (path)
            {
                case 0:x_tmp=x+1;break;
                case 1:x_tmp=x-1;break;
                case 2:y_tmp=y+1;break;
                case 3:y_tmp=y-1;break;
                case 4:z_tmp=z+1;break;
                case 5:z_tmp=z-1;break;                            
            }
            if (check_nucleus(x_tmp,y_tmp,z_tmp)){
                k+=1;
                break;
                }
        }
        x=x_tmp;y=y_tmp;z=z_tmp;
        // cout<<'\t'<<x_tmp<<'\t'<<y_tmp<<'\t'<<z_tmp<<'\n';           
        monomer={x,y,z};
        polymer.push_back(monomer);
        
        //update location 和 contact 
        if (i%4==0){
            if (one_tad_locations[thread_num].find(monomer)!=one_tad_locations[thread_num].end()){
                //如果monomer在location的key中,说明monomer不是第一次出现，说明有两个以上monomer出现contact，此时初始化contact
                for (auto elem : one_tad_locations[thread_num].find(monomer)->second){
                    one_tad_contact[thread_num][{min(elem,position2site_index[i]),max(elem,position2site_index[i])}]=0;//这个0是用来记录后面第几次move  
                }
                one_tad_locations[thread_num].find(monomer)->second.push_back(position2site_index[i]);//三维位置重合的点加入Location
            }
            else { one_tad_locations[thread_num][monomer] = {position2site_index[i]}; }//如果monomer不在location中，即第一次出现，加入该monomer

            //intial polymer_bin
            polymer_bin[thread_num].push_back(monomer);
        }
    }  
}

void initial_polymer_byfile(vector<vector<int>> &polymer,vector<int> &position2site_index,vector<int> &site_index2position,vector<int> &site_available,int pol_length,int thread_num,string file){
    vector<int> monomer;
    int x,y,z,loc;
    int i=0;
    // vector<vector<int>>().swap(polymer);
    site_available.resize(pol_length);
    position2site_index.resize(pol_length);
    site_index2position.resize(pol_length);
    fill(site_available.begin(), site_available.end(), 0);
    fill(position2site_index.begin(), position2site_index.end(), 0);
    fill(site_index2position.begin(), site_index2position.end(), 0);

    ifstream initial_polymer;
    initial_polymer.open(file); 
    
    while(initial_polymer>> x >> y >> z >> loc){

        site_available[loc]=1;
        position2site_index[i]=loc;
        monomer={x,y,z};
        polymer.push_back(monomer); 

        //update location and contact 
        if (loc%6==0){
            if (one_tad_locations[thread_num].find(monomer)!=one_tad_locations[thread_num].end()){
                //如果monomer在location的key中,说明monomer不是第一次出现，说明有两个以上monomer出现contact，此时初始化contact
                for (auto elem : one_tad_locations[thread_num].find(monomer)->second){
                    one_tad_contact[thread_num][{min(elem,position2site_index[i]),max(elem,position2site_index[i])}]=0;//这个0是用来记录后面第几次move  
                }
                one_tad_locations[thread_num].find(monomer)->second.push_back(position2site_index[i]);//三维位置重合的点加入Location
            }
            else { one_tad_locations[thread_num][monomer] = {position2site_index[i]}; }//如果monomer不在location中，即第一次出现，加入该monomer        
        }
        i++;
    } 
    initial_polymer.close();   
    //initial site_index2position
    site_index2position[0]=0;
    for (int i=1;i<pol_length;i++){
        site_index2position[i]=site_available[i]+site_index2position[i-1];
    }    

    vector<vector<int>>().swap(polymer_bin[thread_num]);
    //initial polymer bin
    for(int i=0;i<tad_bin_num*6;i+=6){
        monomer=polymer[site_index2position[i]];        
        polymer_bin[thread_num].push_back(monomer);   
    }
}




//-----------------------tad,refer contact 和tad contact处理--------------------------

//get tad data
void get_tadset(string filepath, vector<vector<int>> &tads_bound){
    ifstream path;
    path.open(filepath);
    int sitei;
    int sitej;
    int index=0;
       
    while(path >> sitei >> sitej ){     
        tads_bound.push_back({sitei,sitej});
        index+=1;
    }          
    path.close();       
}

void initial_oneTad_contact(int bin_start, int bin_end){
    reference_one_tad_contact_s=reference_contact_s.block(bin_start,bin_start,bin_end-bin_start,bin_end-bin_start); 

    Eigen::MatrixXf reference_one_tad_contact=MatrixXf(reference_one_tad_contact_s);

    for (int i=0;i<reference_one_tad_contact.cols();++i){
        reference_one_tad_contact(i,i)=0;
        if(i-1>=0){
            reference_one_tad_contact(i,i-1)=0;
            reference_one_tad_contact(i-1,i)=0;
            }
    }
    reference_one_tad_contact_s=reference_one_tad_contact.sparseView();   
}

void initial_mask(int start_bin,int end_bin){
    mask_bin_list.clear();
    unordered_map<string,int> chr_start_bin_index;
    unordered_map<string,int> chr_end_bin_index;
    string size_file="/data/software/HiC/DATA/reference/hg19/hg19_chr/hg19.chrom_20k.sizes";
    ifstream size;
    size.open(size_file);

    //initial chr_start_bin
    string key;int value;int index=0;
    while(size >> key >> value ){
        chr_start_bin_index[key]=index;
        index+=value;
        chr_end_bin_index[key]=index;
    }   
    
    size.close();  
    int chr_start_bin=chr_start_bin_index[chr];
    int chr_end_bin=chr_end_bin_index[chr];

    //initial mask file
    string mask_file="/data/home/txguo/data_use/maxEnt/reference/mask/20k_mask.txt";
    ifstream mask_f;
    int i;
    mask_f.open(mask_file);  
    while(mask_f >> i ){ 
        if(i>=chr_start_bin && i<chr_end_bin){
            i=i-chr_start_bin;
            if(i>=start_bin && i<end_bin)
                mask_bin_list.push_back(i-start_bin);          
        }
    } 
    mask_f.close();    


    //update mask
    Eigen::MatrixXf mask=Eigen::MatrixXf::Ones(tad_bin_num,tad_bin_num);
    Eigen::MatrixXf zero_row= Eigen::MatrixXf::Zero(1,tad_bin_num);
    Eigen::MatrixXf zero_col= Eigen::MatrixXf::Zero(tad_bin_num,1);
    //mask_s
    for(int elem:mask_bin_list){
        mask.row(elem)=zero_row;
        mask.col(elem)=zero_col;
    }    

    mask_s=mask.sparseView();
    mask_s.selfadjointView<Eigen::Upper>();
    cout<<mask_s.cols();

    reference_one_tad_contact_s=reference_one_tad_contact_s.cwiseProduct(mask_s);
    
    for (int k=0; k<reference_one_tad_contact_s.outerSize(); ++k){
        for (SparseMatrix<float>::InnerIterator it(reference_one_tad_contact_s,k); it; ++it)
        {
            refer_hash[{it.row(),it.col()}]=it.value();
        }
    }      
}

// bool get_bool(int i,int j,vector<int> &index_list){
//     bool res=true;
//     for(auto elem:index_list){
//         res=res & (i!=elem) & (j!=elem);
//     }
//     return res;
// }

// void get_mask(SparseMatrix<float> &matrix){
//     //update mask
//     matrix.prune([&](int i, int j, float) { return get_bool(i,j,mask_bin_list); }); 

//     // for(auto elem:mask_bin_list){
//     //     matrix.row(elem);
//     // }
// }


// void initial_refer(string chr){

//     //initial chrom_size
//     unordered_map<string,int> chr_start_bin_index;
//     string size_file="/data/software/HiC/DATA/reference/hg19/hg19_chr/hg19.chrom_20k.sizes";
//     ifstream size;
//     size.open(size_file);

//     //initial start_bin
//     string key;int value;int index=0;
//     while(size >> key >> value ){
//         chr_index[key]=value;
//         chr_start_bin_index[key]=index;
//         index+=value;
//     }   
    
//     size.close();  


//     // initial the refer data
//     chr_len=chr_index[chr];
//     int start_bin=chr_start_bin_index[chr];
//     MatrixXf reference_contact=MatrixXf::Zero(chr_len,chr_len);

//     // string refer_file="/data/home/txguo/data_use/maxEnt/input_data/normalize/KR_norm/chr"+(chr)+"_20k_corrected_contact_KR.matrix";

//     string refer_file="/data/home/txguo/data_use/data_final/normalize/ICENorm/chr1_icenorm.mtx";
//     // string  refer_file="/data/home/txguo/data_use/maxEnt/reference/reference_use/chr"+chr+"_20k.mtx";
//     loadMarket(reference_contact_s,refer_file);

// }
void initial_refer(string chr){

    //initial chrom_size
    unordered_map<string,int> chr_index;
    unordered_map<string,int> chr_start_bin_index;
    string size_file="/data/software/HiC/DATA/reference/hg19/hg19_chr/hg19.chrom_20k.sizes";
    // string size_file="/data/home/txguo/data_use/maxEnt/reference/KBM7/window/GRCh38.chrom_20k.sizes";
    ifstream size;
    size.open(size_file);

    //initial start_bin
    string key;int value;int index=0;
    while(size >> key >> value ){
        chr_index[key]=value;
        chr_start_bin_index[key]=index;
        index+=value;
    }   
    
    size.close();  


    // initial the refer data
    chr_len=chr_index[chr];
    int start_bin=chr_start_bin_index[chr];
    MatrixXf reference_contact=MatrixXf::Zero(chr_len,chr_len);

    string refer_file="/data/home/txguo/data_use/maxEnt/input_data/normalize/KR_norm/chr"+(chr)+"_20k_corrected_contact_KR.matrix";
    // string refer_file="/data/home/txguo/data_use/maxEnt/reference/KBM7/chr_refer/KBM7_chr"+(chr)+"_20k.mtx";
    ifstream refer;
    refer.open(refer_file);

    int sitei;
    int sitej;
    float contact_num;

    while(refer >> sitei >> sitej >> contact_num){
        if(sitej>sitei+mask_bins){
            reference_contact(sitei,sitej)=contact_num;
        }
    }   
    
    refer.close();  
    reference_contact_s=reference_contact.sparseView();
}

// void initial_refer(string chr){

//     //initial chrom_size
//     unordered_map<string,int> chr_start_bin_index;
//     // string size_file="/data/software/HiC/DATA/reference/hg19/hg19_chr/hg19.chrom_20k.sizes";
//     string size_file="/data/home/txguo/data_use/maxEnt/reference/KBM7/window/GRCh38.chrom_20k.sizes";
//     ifstream size;
//     size.open(size_file);

//     //initial start_bin
//     string key;int value;int index=0;
//     while(size >> key >> value ){
//         chr_index[key]=value;
//         chr_start_bin_index[key]=index;
//         index+=value;
//     }   
    
//     size.close();      

//     // initial the refer data
//     chr_len=chr_index[chr];
//     int start_bin=chr_start_bin_index[chr];
//     cout<<"test:"<<chr_len<<'\n';

//     // string refer_file="/data/home/txguo/data_use/maxEnt/input_data/normalize/KR_norm/chr"+(chr)+"_20k_corrected_contact_KR.matrix";

//     // string refer_file="/data/home/txguo/data_use/data_final/normalize/ICENorm/chr1_icenorm.mtx";
//     // string  refer_file="/data/home/txguo/data_use/maxEnt/reference/reference_use/chr"+chr+"_20k.mtx";
//     string refer_file="/data/home/txguo/data_use/maxEnt/reference/KBM7/chr_refer/KBM7_chr"+(chr)+"_20k.matrix";
//     loadMarket(reference_contact_s,refer_file);
//     MatrixXf reference_contact=MatrixXf(reference_contact_s);

// }

// void initial_tad_energy_byfile(string filepath,int start_bin,int end_bin){
    
//     int bin_length=end_bin-start_bin+1;
//     Eigen::MatrixXf E=Eigen::MatrixXf::Zero(bin_length,bin_length);

//     ifstream energy_file;
//     energy_file.open(filepath);
//     for (int i=0;i<bin_length;i++){
//         for (int j=0;j<bin_length;j++){
//             energy_file >> E(i,j);           
//         }
        
//     }
//     // cout<<E<<'\n';
//     energy_file.close();   

//     Interaction_E_tad.block(start_bin,start_bin,bin_length,bin_length)=E;
// }

void initial_contact_mat(int threshold){
   
    //提前给稀疏矩阵中距离近的占位为0
    tad_final_contacts_s.resize(tad_bin_num,tad_bin_num);
    tad_final_contacts_s.reserve(VectorXi::Constant(tad_bin_num,threshold));
    for(int l=0;l<number_of_threads;l++){
        tad_total_contacts_s[l].resize(tad_bin_num,tad_bin_num);
        tad_total_contacts_s[l].reserve(VectorXi::Constant(tad_bin_num,threshold));
    }
    return ;    
}



void initial_tad_energy(){

    Interaction_E_tad_s.resize(tad_bin_num,tad_bin_num);
    Interaction_E_tad_s.setZero();   

    Interaction_E_tad=Eigen::MatrixXf::Zero(tad_bin_num,tad_bin_num);  

    // for (int i=0;i<tad_bin_num;i++){
    //     for(int k=i+mask_bins;k<tad_bin_num;k++){
    //         if (reference_one_tad_contact(i,k)==0  ) {//no contact
    //             Interaction_E_tad(i,k)=1.0;//the corresponding Eij is set to a high value at the start of the simulation, typically 10
    //             Interaction_E_tad(k,i)=1.0;
    //         }   
    //     }
    // }
}

//----------------------- neighbor deal -----------------------------
void find_neighbor(int x,int y,int z,int site,int m,vector<vector<int>> &polymer,vector<vector<int>> &polymer_bin,vector<int> &site_index2position,vector<int> &tad_neighbor_contact_new,unordered_map<pair<int, int>, float, pair_hash> &neighbor_contact_dis){
    
    // //find neighbor around (x,y,z)

    if (site%6!=0){return;}
    int new_x;int new_y;int new_z;
    int x_up;int y_up;int z_up;int x_down;int y_down;int z_down;float dis;
    vector<int>().swap(tad_neighbor_contact_new);
    unordered_map<pair<int, int>, float, pair_hash>().swap(neighbor_contact_dis);
    
    x_up=x+neighbor_dis;y_up=y+neighbor_dis;z_up=z+neighbor_dis;
    x_down=x-neighbor_dis;y_down=y-neighbor_dis;z_down=z-neighbor_dis;  

    //find neighbor behind site   
    for (int i=0;i<site;i+=6){
        int site_index=i/6;     
        new_x=polymer_bin[site_index][0];new_y=polymer_bin[site_index][1];new_z=polymer_bin[site_index][2];
        if (new_x<=x_up && new_x>=x_down && new_y<=y_up && new_y>=y_down && new_z<=z_up && new_z>=z_down ){
            dis=(sqrt(pow(new_x-x,2)+pow(new_y-y,2)+pow(new_z-z,2)));
            if (dis!=0 && dis<=neighbor_dis){
                tad_neighbor_contact_new.emplace_back(i);
                neighbor_contact_dis[{i,site}]=dis;
            }            
        }
    }

    //find neighbor before site
    for (int i=site+6;i<pol_length-1;i+=6){
        int site_index=i/6;  
        new_x=polymer_bin[site_index][0];new_y=polymer_bin[site_index][1];new_z=polymer_bin[site_index][2];
        if (new_x<=x_up && new_x>=x_down && new_y<=y_up && new_y>=y_down && new_z<=z_up && new_z>=z_down ){
            dis=(sqrt(pow(new_x-x,2)+pow(new_y-y,2)+pow(new_z-z,2)));
            if (dis!=0 && dis<=neighbor_dis){
                tad_neighbor_contact_new.emplace_back(i);
                neighbor_contact_dis[{site,i}]=dis;
            }   
        }     
    }   

    return ;    
}


void find_neighbor_initial(int x,int y,int z,int site,int m,vector<vector<int>> &polymer,vector<vector<int>> &polymer_bin,vector<int> &site_index2position,vector<int> &tad_neighbor_contact_new,unordered_map<pair<int, int>, float, pair_hash> &neighbor_contact_dis){
    
    // //find neighbor around (x,y,z)

    if (site%6!=0){return;}
    int new_x;int new_y;int new_z;
    int x_up;int y_up;int z_up;int x_down;int y_down;int z_down;float dis;
    vector<int>().swap(tad_neighbor_contact_new);
    unordered_map<pair<int, int>, float, pair_hash>().swap(neighbor_contact_dis);
    
    x_up=x+neighbor_dis;y_up=y+neighbor_dis;z_up=z+neighbor_dis;
    x_down=x-neighbor_dis;y_down=y-neighbor_dis;z_down=z-neighbor_dis;  

    //find neighbor behind site   
    for (int i=0;i<site;i+=6){
        int site_index=i/6;     
        new_x=polymer_bin[site_index][0];new_y=polymer_bin[site_index][1];new_z=polymer_bin[site_index][2];
        dis=(sqrt(pow(new_x-x,2)+pow(new_y-y,2)+pow(new_z-z,2)));
        if (dis!=0 && dis<=20){
            tad_neighbor_contact_new.emplace_back(i);
            neighbor_contact_dis[{i,site}]=dis;
        }            
    }

    //find neighbor before site
    for (int i=site+6;i<pol_length-1;i+=6){
        int site_index=i/6;  
        new_x=polymer_bin[site_index][0];new_y=polymer_bin[site_index][1];new_z=polymer_bin[site_index][2];
        dis=(sqrt(pow(new_x-x,2)+pow(new_y-y,2)+pow(new_z-z,2)));
        if (dis!=0 && dis<=20){
            tad_neighbor_contact_new.emplace_back(i);
            neighbor_contact_dis[{site,i}]=dis;
        }   
    }   

    return ;    
}
// void find_neighbor(int x,int y,int z,int site,int m,vector<vector<int>> &polymer,vector<vector<int>> &polymer_bin,vector<int> &site_index2position,vector<int> &tad_neighbor_contact_new,unordered_map<pair<int, int>, float, pair_hash> &neighbor_contact_dis){
    
//     //find neighbor around (x,y,z)

//     if (site%6!=0){return;}
//     int new_x;int new_y;int new_z;
//     int x_up;int y_up;int z_up;int x_down;int y_down;int z_down;float dis;
//     vector<int>().swap(tad_neighbor_contact_new);
//     unordered_map<pair<int, int>, float, pair_hash>().swap(neighbor_contact_dis);
    
//     x_up=x+neighbor_dis;y_up=y+neighbor_dis;z_up=z+neighbor_dis;
//     x_down=x-neighbor_dis;y_down=y-neighbor_dis;z_down=z-neighbor_dis;  

//     //find neighbor behind site   
//     for (int i=0;i<site;i+=6){
//         int site_index=i/6;     
//         new_x=polymer_bin[site_index][0];new_y=polymer_bin[site_index][1];new_z=polymer_bin[site_index][2];
//         if (new_x<=x_up && new_x>=x_down && new_y<=y_up && new_y>=y_down && new_z<=z_up && new_z>=z_down ){
//             dis=(sqrt(pow(new_x-x,2)+pow(new_y-y,2)+pow(new_z-z,2)));
//             if (dis!=0 && dis<=neighbor_dis){
//                 tad_neighbor_contact_new.emplace_back(i);
//                 neighbor_contact_dis[{i,site}]=dis;
//             }            
//         }
//     }

//     //find neighbor before site
//     for (int i=site+6;i<pol_length-1;i+=6){
//         int site_index=i/6;  
//         new_x=polymer_bin[site_index][0];new_y=polymer_bin[site_index][1];new_z=polymer_bin[site_index][2];
//         if (new_x<=x_up && new_x>=x_down && new_y<=y_up && new_y>=y_down && new_z<=z_up && new_z>=z_down ){
//             dis=(sqrt(pow(new_x-x,2)+pow(new_y-y,2)+pow(new_z-z,2)));
//             if (dis!=0 && dis<=neighbor_dis){
//                 tad_neighbor_contact_new.emplace_back(i);
//                 neighbor_contact_dis[{site,i}]=dis;
//             }   
//         }     
//     }   

//     return ;    
// }


void initial_neighbor(int thread_num){

    //initial tad_neighbor_contact
    // vector<unordered_map<int, vector<int>>>(number_of_threads).swap(tad_neighbor_contact); 
    unordered_map<int, vector<int>>().swap(tad_neighbor_contact[thread_num]); 
    unordered_map<pair<int, int>, float, pair_hash>().swap(neighbor_contact_dis[thread_num]);

    int x;int y;int z;
    vector<int> tad_neighbor_contact_initial;       

    for (int i=0;i<polymer_bin[thread_num].size();i++){
        x=polymer_bin[thread_num][i][0];y=polymer_bin[thread_num][i][1];z=polymer_bin[thread_num][i][2];
        // if(thread_num==0){
        //     cout<<x<<'\t'<<y<<'\t'<<z<<'\n';
        // }
        
        find_neighbor_initial(x,y,z,6*i,0,polymer[thread_num],polymer_bin[thread_num],site_index2position[thread_num],tad_neighbor_contact_initial,neighbor_contact_dis[thread_num]);
 
        //update
        if (tad_neighbor_contact_initial.empty()==0){
            tad_neighbor_contact[thread_num][i*6]=tad_neighbor_contact_initial;

        }    
        for(auto elem:neighbor_contact_dis[thread_num]){
            int sitei=elem.first.first;
            int sitej=elem.first.second;
            float dis=elem.second;
            
            neighbor_contact_data[thread_num][{min(sitei,sitej),max(sitei,sitej)}].dis=dis;
        }        
    }

    if(thread_num==0){
        cout<<"test"<<neighbor_contact_data[thread_num].size()<<'\n';
    }
    
}        


#endif