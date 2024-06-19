#ifndef INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H
#define INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "global.h"
using namespace std;
using namespace Eigen;

//============================= tads calculate ============================

//----------------------------tad_level------------------------
float get_disAenergy(int thread_num,int chr_homology,int tad_index_start,int tad_index,int chr_neighbor_homology,int dir,vector<float> &new_center_vec,vector<MatrixXf> &tads_center_chr,vector<vector<float>> &homology_dis_matrix){

    //dir=1 tad_index>tad_mid;  dir=0 tad_index<tad_mid  
    int tad_len=chr_tad_len[chr_homology/2];
    int change_index_len=homology_dis_matrix_change[thread_num].size();
    int tad_neighbor_chr_len=chr_tad_len[chr_neighbor_homology/2];
    int indexi,indexj,indexi_haploid,indexj_haploid;
    float x,y,z,dis,dis_pre,energy_change=0;
    float energy_change1=0,energy_change2=0;

    int homology_start = homology_tad_start[chr_homology];
    int neighbor_homology_start = homology_tad_start[chr_neighbor_homology];   

    int haploid_start = haploid_tad_start[chr_homology/2];
    int neighbor_haploid_start = haploid_tad_start[chr_neighbor_homology/2]; 

    //change neighbor_chr_tad_center 2 tad_center_vec
    vector<float> tad_center_vec(tads_center_chr[chr_neighbor_homology].data(),tads_center_chr[chr_neighbor_homology].data()+tads_center_chr[chr_neighbor_homology].size());
    for(int tad_neighbor_index=0;tad_neighbor_index<tad_neighbor_chr_len;tad_neighbor_index++){
        
        if(reconstruct_chrom_mark==1 && tad_neighbor_index!=chr_homology){continue;}
        //update distance
        if(dir==1){
            indexi=tad_index-tad_index_start;            
            if(chr_neighbor_homology==chr_homology && tad_neighbor_index>=tad_index_start ){
                x=new_center_vec[change_index_len*0+tad_neighbor_index-tad_index_start]-new_center_vec[change_index_len*0+indexi];   
                y=new_center_vec[change_index_len*1+tad_neighbor_index-tad_index_start]-new_center_vec[change_index_len*1+indexi];
                z=new_center_vec[change_index_len*2+tad_neighbor_index-tad_index_start]-new_center_vec[change_index_len*2+indexi];

            }    
            else{
                x=tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index]-new_center_vec[change_index_len*0+indexi];   
                y=tad_center_vec[tad_neighbor_chr_len*1+tad_neighbor_index]-new_center_vec[change_index_len*1+indexi];
                z=tad_center_vec[tad_neighbor_chr_len*2+tad_neighbor_index]-new_center_vec[change_index_len*2+indexi];                 
            } 
        }
        else if (dir==0){           
            indexi=tad_index;   
            if(chr_neighbor_homology==chr_homology && tad_neighbor_index<=tad_index_start ){
                x=new_center_vec[change_index_len*0+tad_neighbor_index]-new_center_vec[change_index_len*0+indexi];   
                y=new_center_vec[change_index_len*1+tad_neighbor_index]-new_center_vec[change_index_len*1+indexi];
                z=new_center_vec[change_index_len*2+tad_neighbor_index]-new_center_vec[change_index_len*2+indexi]; 
            }    
            else{
                x=tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index]-new_center_vec[change_index_len*0+indexi];   
                y=tad_center_vec[tad_neighbor_chr_len*1+tad_neighbor_index]-new_center_vec[change_index_len*1+indexi];
                z=tad_center_vec[tad_neighbor_chr_len*2+tad_neighbor_index]-new_center_vec[change_index_len*2+indexi];                 
            }             
        }       

        dis=sqrt(x*x+y*y+z*z);


        indexj=neighbor_homology_start+tad_neighbor_index;
        
        dis_pre=homology_dis_matrix[min(homology_start+tad_index,indexj)][max(homology_start+tad_index,indexj)];   
        
        homology_dis_matrix_change[thread_num][min(indexi,indexj)][max(indexi,indexj)]=dis;             

        //update energy   
        indexi_haploid=haploid_start+tad_index;
        indexj_haploid=neighbor_haploid_start+tad_neighbor_index;  
        float pre_e=energy_change;

        if(mark==0){
            if(chr_homology/2==chr_neighbor_homology/2){
                energy_change -= 100*chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis_pre,float(0.5));  
                energy_change += 100*chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis,float(0.5));  
                // energy_change -= chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis_pre,float(0.5));  
                // energy_change += chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis,float(0.5));         
            }
            else{
                energy_change -= chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis_pre,float(1))/max(dis_pre,float(1));  
                energy_change += chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis,float(1))/max(dis,float(1));              
            }            
        }
        else{
            if(chr_homology/2==chr_neighbor_homology/2){
                // energy_change -= 100*chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis_pre,float(0.5));  
                // energy_change += 100*chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis,float(0.5));  
                energy_change -= chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis_pre,float(0.5));  
                energy_change += chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis,float(0.5));         
            }
            else{
                energy_change -= 5*chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis_pre,float(1))/max(dis_pre,float(1));  
                energy_change += 5*chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis,float(1))/max(dis,float(1));              
            }             
        }


 
    } 

    return energy_change;              

}


float get_disAenergy2chrom(int thread_num,int chr_homology,int tad_index_start,int tad_index,int chr_neighbor_homology,int dir,vector<float> &new_center_vec,vector<MatrixXf> &tads_center_chr,vector<vector<float>> &homology_dis_matrix){

    //dir=1 tad_index>tad_mid;  dir=0 tad_index<tad_mid  
    int change_index_len;
    int tad_len=chr_tad_len[chr_homology/2];
    int tad_neighbor_chr_len=chr_tad_len[chr_neighbor_homology/2];
    int indexi,indexj,indexi_haploid,indexj_haploid;
    float x=0,y=0,z=0,dis,dis_pre,energy_change=0;

    int homology_start = homology_tad_start[chr_homology];
    int neighbor_homology_start = homology_tad_start[chr_neighbor_homology];   

    int haploid_start = haploid_tad_start[chr_homology/2];
    int neighbor_haploid_start = haploid_tad_start[chr_neighbor_homology/2]; 

    // cout<<change_index_len<<'\t'<<tad_index_start<<'\t'<<tad_index<<'\n';

    //change neighbor_chr_tad_center 2 tad_center_vec

    vector<float> tad_center_vec(tads_center_chr[chr_neighbor_homology].data(),tads_center_chr[chr_neighbor_homology].data()+tads_center_chr[chr_neighbor_homology].size());
    for(int tad_neighbor_index=0;tad_neighbor_index<tad_neighbor_chr_len;tad_neighbor_index++){

        //update distance   
        indexi=tad_index;
        indexj=tad_neighbor_index;           
    
        if(dir==1){
            change_index_len=tad_len-tad_index_start;
            if(chr_neighbor_homology==chr_homology && tad_neighbor_index>=tad_index_start ){
                x=new_center_vec[change_index_len*0+tad_neighbor_index-tad_index_start]-new_center_vec[change_index_len*0+tad_index-tad_index_start];   
                y=new_center_vec[change_index_len*1+tad_neighbor_index-tad_index_start]-new_center_vec[change_index_len*1+tad_index-tad_index_start];
                z=new_center_vec[change_index_len*2+tad_neighbor_index-tad_index_start]-new_center_vec[change_index_len*2+tad_index-tad_index_start];
            }    
            else{
                x=tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index]-new_center_vec[change_index_len*0+tad_index-tad_index_start];   
                y=tad_center_vec[tad_neighbor_chr_len*1+tad_neighbor_index]-new_center_vec[change_index_len*1+tad_index-tad_index_start];
                z=tad_center_vec[tad_neighbor_chr_len*2+tad_neighbor_index]-new_center_vec[change_index_len*2+tad_index-tad_index_start];                 
            }                  
        }
        else if (dir==0){ 
            change_index_len=tad_index_start+1;
            if(chr_neighbor_homology==chr_homology && tad_neighbor_index<=tad_index_start ){
                x=new_center_vec[change_index_len*0+tad_neighbor_index]-new_center_vec[change_index_len*0+indexi];   
                y=new_center_vec[change_index_len*1+tad_neighbor_index]-new_center_vec[change_index_len*1+indexi];
                z=new_center_vec[change_index_len*2+tad_neighbor_index]-new_center_vec[change_index_len*2+indexi]; 
            }    
            else{
                x=tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index]-new_center_vec[change_index_len*0+indexi];   
                y=tad_center_vec[tad_neighbor_chr_len*1+tad_neighbor_index]-new_center_vec[change_index_len*1+indexi];
                z=tad_center_vec[tad_neighbor_chr_len*2+tad_neighbor_index]-new_center_vec[change_index_len*2+indexi];                 
            }   
        }

        // if(isnan(x))cout<<"x:"<<x<<'\t'<<indexi<<'\t'<<new_center_vec.size()<<'\t'<<new_center_vec[change_index_len*0+indexi]<<'\n';          
             
        dis=sqrt(x*x+y*y+z*z);      
        // if(dis>50 && abs(tad_index-tad_neighbor_index)<3){cout<<"dis:"<<dis<<'\t'<<tad_index<<'\t'<<tad_neighbor_index<<'\t'<<x<<'\t'<<y<<'\t'<<z<<'\t'<<tad_center_vec[tad_neighbor_chr_len*2+tad_neighbor_index]<<'\t'<<new_center_vec[change_index_len*2+tad_index-tad_index_start]<<'\n';}
        dis_pre=homology_dis_matrix[min(indexi,indexj)][max(indexi,indexj)];   

        homology_dis_matrix_change[thread_num][min(indexi,indexj)][max(indexi,indexj)]=dis;             

        //update energy   

        if(chr_homology/2==chr_neighbor_homology/2){
            energy_change -= 100*chr_E[min(indexi,indexj)][max(indexi,indexj)]/pow(max(dis_pre,float(1)),1);  
            energy_change += 100*chr_E[min(indexi,indexj)][max(indexi,indexj)]/pow(max(dis,float(1)),1);  
        }
        else{
            energy_change -= chr_E[min(indexi,indexj)][max(indexi,indexj)]/max(dis_pre,float(1))/max(dis_pre,float(1));  
            energy_change += chr_E[min(indexi,indexj)][max(indexi,indexj)]/max(dis,float(1))/max(dis,float(1));              
        }   

        // if(energy_change!=0)cout<<"dis:"<<dis<<'\t'<<dis_pre<<'\t'<<energy_change<<'\t'<<indexi<<'\t'<<indexj<<'\n';  
  
        
    } 
    return energy_change;              

}

float get_delta_E(int thread_num,int chr_homology,int tad_index,MatrixXf tads_center_change,vector<MatrixXf> &tads_center_chr,vector<vector<float>> &homology_dis_matrix){

    int chr_another_homology,another_homology_index,tad_neighbor_chr_len,indexi,indexj,dis;
    int tad_len=chr_tad_len[chr_homology/2];
    int tad_start=homology_tad_start[chr_homology];
    float energy_change=0;

    MatrixXf new_center;

    VectorXi center,point,margin;
    vector<float> new_center_vec(tads_center_change.data(),tads_center_change.data()+tads_center_change.size());    
    int change_index_len=tad_len-tad_index;
    // cout<<"x:"<<tad_index<<'\t'<<new_center_vec[change_index_len*0+tad_index-tad_index]<<'\t'<<tads_center_chr[chr_homology](tad_index,0)<<'\n';
    // cout<<"y:"<<new_center_vec[change_index_len*1+tad_index-tad_index]<<'\t'<<tads_center_chr[chr_homology](tad_index,1)<<'\n';
    // cout<<"z:"<<new_center_vec[change_index_len*2+tad_index-tad_index]<<'\t'<<tads_center_chr[chr_homology](tad_index,2)<<'\n';    
    if(tad_index>=tad_len/2){

        for(int i=tad_index;i<tad_len;i++){    
            for(int chr_neighbor=0;chr_neighbor<chr_homology_num;chr_neighbor++){ 
                if(reconstruct_chrom_mark==1 && chr_neighbor!=chr_homology){continue;} 
                energy_change+=get_disAenergy2chrom(thread_num,chr_homology,tad_index,i,chr_neighbor,1,new_center_vec,tads_center_chr,homology_dis_matrix); 
                // if(energy_change!=0){cout<<"energy:"<<energy_change<<'\n';}
            }
        }     
    }
    else if(tad_index<tad_len/2 ){
    
        for(int i=0;i<tad_index+1;i++){    
            for(int chr_neighbor=0;chr_neighbor<chr_homology_num;chr_neighbor++){ 
                if(reconstruct_chrom_mark==1 && chr_neighbor!=chr_homology){continue;} 
                energy_change+=get_disAenergy2chrom(thread_num,chr_homology,tad_index,i,chr_neighbor,0,new_center_vec,tads_center_chr,homology_dis_matrix); 
                // if(energy_change!=0){cout<<"energy:"<<energy_change<<'\n';}
            }
        }
    } 
    

    // if(energy_change!=0){cout<<"energy:"<<energy_change<<'\n';}

    return energy_change;
}

void update_dis_matrix(int thread_num,int chr_homology,int tad_index,int dir){

    //update neighbor
    int indexi,indexj;
    int chr_len=chr_tad_len[chr_homology/2];

    for(int i=0;i<chr_len;i++){
        for(int j=i+1;j<chr_len;j++){              
            homology_dis_matrix[thread_num][i][j]=homology_dis_matrix_change[thread_num][i][j];
        }
    }
     


}


// void update_dis_matrix(int thread_num,int chr_homology,int tad_index,int dir){

//     //update neighbor
//     int indexi,indexj;
//     if(dir==1){//tad_index>=mid_tad_index
//         vector<vector<float>> change_dis;
//         int change_index=homology_tad_start[chr_homology]+tad_index;
//         int end_index=homology_tad_start[chr_homology]+chr_tad_len[chr_homology/2];
        
//         for(int i=0;i<homology_dis_matrix_change[thread_num].size();i++){
//             for(int j=0;j<tads_num;j++){
//                 indexi=change_index+i;
//                 indexj=j;
                  
//                 if(indexi>indexj){
//                     swap(indexi,indexj);
//                 }
             
//                 homology_dis_matrix[thread_num][indexi][indexj]=homology_dis_matrix_change[thread_num][i][j];
//             }
//         }
//     }
//     else if(dir==0){//tad_index<mid_tad_index 
//         vector<vector<float>> change_dis;
//         int change_index=homology_tad_start[chr_homology];

//         for(int i=0;i<homology_dis_matrix_change[thread_num].size();i++){
//             for(int j=0;j<tads_num;j++){
//                 indexi=change_index+i;
//                 indexj=j;
//                 if(indexi>indexj){
//                     swap(indexi,indexj);            
//                 }
//                 homology_dis_matrix[thread_num][indexi][indexj]=homology_dis_matrix_change[thread_num][i][j];
//             }
//         }        
//     }
// }

void update_neighbor_contact(int homology,int thread_num,int m){
    int m_pre=homology_m[thread_num][homology];
    int chr_len=chr_tad_len[homology/2];
    int homology_start=homology_tad_start[homology];
    int homology_end=homology_start+chr_len;
    float dis;
    float contact;


    // m=thread_m[thread_num];
    if(m==m_pre){return;}
    
    for(int i=0;i<chr_len;i++){
        for(int j=i+1;j<chr_len;j++){
            // if(j>=homology_start && j<homology_end){continue;}
            dis=max(homology_dis_matrix[thread_num][i][j],float(1));
            // contact=float(m-m_pre)/pow(dis,2.8);
            contact=float(m-m_pre)/pow(dis,2.8);
            homology_contact_matrix[thread_num][i][j]+=contact;
        }
    }
 

    // for(int i=homology_start;i<homology_end;i++){
    //     for(int j=0;j<tads_num;j++){
    //         // if(j>=homology_start && j<homology_end){continue;}
    //         dis=max(homology_dis_matrix[thread_num][min(i,j)][max(i,j)],float(1));
    //         contact=float(m-m_pre)/pow(dis,2.4);
    //         homology_contact_matrix[thread_num][min(i,j)][max(i,j)]+=contact;

    //     }    

    //     // //防止染色体内contact自加两次
    //     // for(int j=i+1;j<homology_end;j++){
            
    //     //     dis=max(homology_dis_matrix[thread_num][i][j],float(1));
    //     //     contact=float(m-m_pre)/pow(dis,1.8);
    //     //     homology_contact_matrix[thread_num][i][j]+=contact;

    //     // }               
    // }        

    homology_m[thread_num][homology]=m;
}


//------------------------update tad contact-----------------------------
void update_tads_location(MatrixXi tad ,MatrixXi change_tad ,int tad_start_index,int tad_end_index,int thread_num,int m){

    vector<int> monomer;
    vector<int> change_monomer;
    int start;int end;int last;int loc;float dis;
    int first_index=tads_bound[0][0]; 
    int first=tads_bound[tad_start_index][0];
    if (tad_end_index==tads_bound.size()-1){last=tads_bound[tad_end_index][1];}  
    else{last=tads_bound[tad_end_index][1];}
    //Throw away old contacts, at the same time update contact frequency map 更新旧的site位点接触的contact
 
    //update tads location
    for (int site=first;site<last;++site){
        loc=site - first;
        monomer = {tad(loc,0),tad(loc,1),tad(loc,2)};
        change_monomer = {change_tad(loc,0),change_tad(loc,1),change_tad(loc,2)};

        //erase old contacts
        for (auto elem : tads_location[thread_num][monomer]){
            if (elem < first || elem >=last ){     
                tads_total_contacts_s[thread_num].coeffRef(tad_bin_dict[min(elem,site)],tad_bin_dict[max(elem,site)]) += m - tads_contact[thread_num][{min(elem,site),max(elem,site)}];
                tads_contact[thread_num].erase({min(elem,site),max(elem,site)});
            }
        } 

        //put in new contacts
        if (tads_location[thread_num].find(change_monomer)!=tads_location[thread_num].end()){//新site与location中位置重合
            for (auto elem : tads_location[thread_num][change_monomer]) {
                if (elem < first || elem >=last){                
                    tads_contact[thread_num][{min(elem,site),max(elem,site)}] = m;
                }
            }
        }  

        //=========update hash map locations=========
        //erase old site
        if (tads_location[thread_num][monomer].size()==1){
            tads_location[thread_num].erase(monomer);
        }
        else {
            tads_location[thread_num][monomer].erase(find(tads_location[thread_num][monomer].begin(), tads_location[thread_num][monomer].end(),site));
        }
        //add new site
        if (tads_location[thread_num].find(change_monomer) != tads_location[thread_num].end()){
            tads_location[thread_num][change_monomer].emplace_back(site);
        }
        else {
            tads_location[thread_num][change_monomer] = {site};
        }                 
    }   
        
}

void update_tads_change_location(MatrixXi tad ,MatrixXi change_tad ,int tad_start_index,int thread_num,int m){

    vector<int> monomer;
    vector<int> change_monomer;
    int start;int end;int loc;float dis;
    int first=tads_bound[tad_start_index][0]; 
    int first_index=tads_bound[0][0];
    int last=tads_bound[tads_bound.size()-1][1];  

    // //update neighbor
    // for(int i=0;i<tads_bound.size();i++){
    //     for(int j=i+1;j<tads_bound.size();j++){
    //         dis=tads_center_dis[thread_num](i,j);
    //         if(i!=j ){
    //             tads_total_contacts[thread_num](min(i,j),max(i,j)) +=min(pow(dis,-3),1.00);        
    //         }
    //     }
    // } 


    //Throw away old contacts, at the same time update contact frequency map 更新旧的site位点接触的contact

    for (int site=first;site<last;++site){
        loc=site - first;
        monomer = {tad(loc,0),tad(loc,1),tad(loc,2)};
        change_monomer = {change_tad(loc,0),change_tad(loc,1),change_tad(loc,2)};
        start=tads_bound[tad_bin_dict[site]][0];
        end=tads_bound[tad_bin_dict[site]][1];

        //erase old contacts
        for (auto elem : tads_location[thread_num][monomer]){
            if (elem < start ||  elem >=end){                      
                tads_total_contacts_s[thread_num].coeffRef(tad_bin_dict[min(elem,site)],tad_bin_dict[max(elem,site)]) += m - tads_contact[thread_num][{min(elem,site),max(elem,site)}]; ;
                tads_contact[thread_num].erase({min(elem,site),max(elem,site)});
            }
        } 

        //put in new contacts
        if (tads_location[thread_num].find(change_monomer)!=tads_location[thread_num].end()){//新site与location中位置重合
            for (auto elem : tads_location[thread_num][change_monomer]) {
                if (elem < start ||  elem >=end){               
                    tads_contact[thread_num][{min(elem,site),max(elem,site)}] = m;
                }
            }
        }  

        //=========update hash map locations=========
        //erase old site
        if (tads_location[thread_num][monomer].size()==1){
            tads_location[thread_num].erase(monomer);
        }
        else {
            tads_location[thread_num][monomer].erase(find(tads_location[thread_num][monomer].begin(), tads_location[thread_num][monomer].end(),site));
        }
        //add new site
        if (tads_location[thread_num].find(change_monomer) != tads_location[thread_num].end()){
            tads_location[thread_num][change_monomer].emplace_back(site);
        }
        else {
            tads_location[thread_num][change_monomer] = {site};
        }                 
    }  
   
}


void update_tads_center(int index,MatrixXi changed_structure,MatrixXf &changed_tads_center){
    MatrixXi x;MatrixXi y;MatrixXi z;float x_c;float y_c;float z_c;int start;int end;
    int first=tads_bound[index][0];
    for (int i=index;i<tads_bound.size();++i){
        start=tads_bound[i][0]-first;
        end=tads_bound[i][1]-first;

        x=changed_structure.middleRows(start,end-start).col(0);
        y=changed_structure.middleRows(start,end-start).col(1);
        z=changed_structure.middleRows(start,end-start).col(2);
        x_c=float(x.sum())/(end-start);
        y_c=float(y.sum())/(end-start);
        z_c=float(z.sum())/(end-start);
        changed_tads_center.row(i-index)<<x_c,y_c,z_c;          
    }
}


void update_tads_center_data(int thread,int index,int m){
    Eigen::MatrixXf margin;
    for (int i=0;i<tads_center[thread].rows();++i){
        margin=tads_center[thread].row(index)-tads_center[thread].row(i);
        tads_center_dis[thread](i,index)=sqrt((margin.cwiseProduct(margin)).sum());    
        tads_center_dis[thread](index,i)=tads_center_dis[thread](i,index);

        tads_center_m[thread](index,i)=m;
        tads_center_m[thread](i,index)=m;
    }
}

void normalize_tads() {   

    float sum ;
    // int length=tads_num/2; 
    int length=reference_genome_tads_contact.size();
    for(int i=0;i<length;i++){

        //mask the map
        if(accumulate(reference_genome_tads_contact[i].begin(),reference_genome_tads_contact[i].end(),float(0))==0){
            for (int j = 0; j < length; j++) {
                haploid_contact_matrix[i][j]=0;
                haploid_contact_matrix[j][i]=0;
            }
            continue;
        }    
          
        //calculate sum  
        for(int j=i+1;j<length;j++){
            sum+=haploid_contact_matrix[i][j];
        }
    }
    float multiply= float(length)/(2*sum);       

    for(int i=0;i<length;i++){
        for(int j=i+1;j<length;j++){
            haploid_contact_matrix[i][j]=haploid_contact_matrix[i][j]*multiply;
        }
    }    

    sum=0;
    for(int i=0;i<length;i++){
        for(int j=i+1;j<length;j++){
            sum+=haploid_contact_matrix[i][j];
        }
    }   

    cout<<"nor_check"<<"\t"<<sum<<"\n";
}

void normalize_tads1() {   

    float sum_intra=0,sum_inter=0,sum_intra_refer=0,sum_inter_refer=0,sum=0,sum_refer=0;
    int length=tads_num/2; 

    for(int i=0;i<length;i++){
        //mask the map
        if(accumulate(reference_genome_tads_contact[i].begin(),reference_genome_tads_contact[i].end(),float(0))==0){
            for (int j = 0; j < length; j++) {
                haploid_contact_matrix[i][j]=0;
                haploid_contact_matrix[j][i]=0;
            }
            continue;
        }    
        //calculate sum        
        for(int j=i+1;j<length;j++){
            sum_refer+=reference_genome_tads_contact[i][j];
            sum+=haploid_contact_matrix[i][j];
        }
    }

    for(int i=0;i<22;i++){
        int chr_start=haploid_tad_start[i];
        int chr_end=chr_tad_len[i]+chr_start;

        for(int j=chr_start;j<chr_end;j++){
            for(int k=j+1;k<chr_end;k++){
                sum_intra+=haploid_contact_matrix[j][k];
                sum_intra_refer+=reference_genome_tads_contact[j][k];
            }
        }
    }

    sum_inter_refer=sum_refer-sum_intra_refer;
    sum_inter=sum-sum_intra;
    float multiple_intra=sum_refer/sum;
    float multiple_inter=sum_inter_refer/sum_inter;
    cout<<"multiple_intra:"<<multiple_intra<<'\n';
    cout<<"multiple_inter:"<<multiple_inter<<'\n';

    //normalize
    for(int i=0;i<length;i++){
        for(int j=i+1;j<length;j++){
            haploid_contact_matrix[i][j]=haploid_contact_matrix[i][j]*multiple_inter;
        }
    }    

    for(int i=0;i<22;i++){
        int chr_start=haploid_tad_start[i];
        int chr_end=chr_tad_len[i]+chr_start;

        for(int j=chr_start;j<chr_end;j++){
            for(int k=j+1;k<chr_end;k++){
                haploid_contact_matrix[j][k]=haploid_contact_matrix[j][k]/multiple_inter*multiple_intra;
            }
        }
    }    

    

    sum=0;
    for(int i=0;i<length;i++){
        for(int j=i+1;j<length;j++){
            sum+=haploid_contact_matrix[i][j];
        }
    }   

    cout<<"nor_check"<<"\t"<<sum<<"\n";
}

void normalize_tads_reference(SparseMatrix<float> &reference_tads_contact_s) {

    int length=reference_tads_contact_s.cols();
    float sum=reference_tads_contact_s.sum();
    reference_tads_contact_s*=float(length)/(2*sum);
    cout<<"refer_check"<<"\t"<<reference_tads_contact_s.sum()<<"\n";

    MatrixXf reference_genome_tads_contact_dense=MatrixXf(reference_tads_contact_s);

    //initial reference_genome_tads_contact
    (vector<vector<float>>(length,vector<float>(length,0))).swap(reference_genome_tads_contact);
    for(int i=0;i<length;i++){
        for(int j=i+1;j<length;j++){
            reference_genome_tads_contact[i][j]=reference_genome_tads_contact_dense(i,j);
            reference_genome_tads_contact[j][i]=reference_genome_tads_contact[i][j];
        }
    }   

}

void update_tads_energies() {
    //update energy & combine real contact and refer contact
    float checksum = 0;
    float loss=0;
    // int length=tads_num/2;
    int length=reference_genome_tads_contact.size();
    cout<<"length"<<'\t'<<length<<'\n';
    float margin;

    for (int i = 0; i < length; i++) {

        //mask the map
        if(accumulate(reference_genome_tads_contact[i].begin(),reference_genome_tads_contact[i].end(),float(0))==0){
            for (int j = 0; j < length; j++) {
                haploid_contact_matrix[i][j]=0;
                haploid_contact_matrix[j][i]=0;
            }
            continue;
        }
        //update energy
        for (int j = i+1; j < length; j++) {
            margin=haploid_contact_matrix[i][j]-reference_genome_tads_contact[i][j];         
            chr_E[i][j] += learning_rate * sqrt(1/(max(reference_genome_tads_contact[i][j],float(0.0001))))*margin;
            chr_E[j][i] = chr_E[i][j];
            loss+=fabs(margin);        

            //combine real contact and refer contact
            haploid_contact_matrix[j][i]=reference_genome_tads_contact[i][j];
        }
    } 
    cout<<"loss"<<"\t"<<loss<<"\n";
}

void get_new_tads_energies() {
    //update energy & combine real contact and refer contact
    float checksum = 0;
    float loss=0;
    // int length=tads_num/2;
    int length=reference_genome_tads_contact.size();
    cout<<"length"<<'\t'<<length<<'\n';
    float margin;

    for (int i = 0; i < length; i++) {

        //mask the map
        if(accumulate(reference_genome_tads_contact[i].begin(),reference_genome_tads_contact[i].end(),float(0))==0){
            for (int j = 0; j < length; j++) {
                haploid_contact_matrix[i][j]=0;
                haploid_contact_matrix[j][i]=0;
            }
            continue;
        }
        //update energy
        for (int j = i+1; j < length; j++) {
            margin=haploid_contact_matrix[i][j]-reference_genome_tads_contact[i][j];         
            chr_E[i][j] = learning_rate * sqrt(1/(max(reference_genome_tads_contact[i][j],float(0.0001))))*margin;
            chr_E[j][i] = chr_E[i][j];
            loss+=fabs(margin);        

            //combine real contact and refer contact
            haploid_contact_matrix[j][i]=reference_genome_tads_contact[i][j];
        }
    } 
    cout<<"loss"<<"\t"<<loss<<"\n";
 
}


void normalize() {   
    int first =tads_bound[0][0];
    int length=final_contacts.rows();
    float sum = final_contacts.sum();
    float check1=0;
    float multiply= (length)/(2*sum);

    final_contacts*= multiply;// pol_length/4     /(8*sum)=/(4*4*sum)*2
    check1+=final_contacts.sum();         
    cout<<"nor_check"<<"\t"<<check1<<"\n";
}

void normalize_reference() {
    int first =tads_bound[0][0];
    int length=refer_contact.rows();
    double sum = refer_contact.sum();
    refer_contact*=(length)/(2*sum);
    
    cout<<"refer_check"<<"\t"<<refer_contact.sum()<<"\n";
}

// void update_energies() {
//     float checksum = 0;
//     float loss=0;
//     int first =tads_bound[0][0];
//     int length=refer_contact.rows();
//     Eigen::MatrixXd margin;
//     Eigen::MatrixXd One=Eigen::MatrixXd::Ones(length,length) ;
//     margin=final_contacts - refer_contact;
//     for (int i = 0; i < length; i++) {
//         for (int j = i+1; j < length; j++) {

//             Interaction_E(i,j) += learning_rate * sqrt(1/(max(refer_contact(i,j),double(0.00001))))*(margin(i,j));
//             // cout<<sqrt(1/(max(refer_contact(i,j),double(0.00001))))*(double(final_contacts(i,j))- refer_contact(i,j));
//             Interaction_E(j,i) = Interaction_E(i,j);
//             loss+=fabs(margin(i,j));
//             checksum += (final_contacts(i,j))*length;//checksum 512 (32*32/2)                        
//         }
//     }
//     cout<<"loss"<<"\t"<<loss<<"\n";

//     //calculates the shift of all energies, imposed to ensure a MaxEnt solution is found for the contact frequency scale
    
//     float shift=(2*(Interaction_E.cwiseProduct(refer_contact)) / length).sum();

//     cout << "Shift: " << shift << endl;
//     //update E
//     Interaction_E-=shift*One;
// }

#endif //INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H