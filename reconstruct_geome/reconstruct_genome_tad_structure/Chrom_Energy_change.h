#ifndef CHROM_ENERGY_CHANGES_H
#define CHROM_ENERGY_CHANGES_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "global.h"
#include "initial_polymer.h"

using namespace std;
using namespace Eigen;


//------------------------------chromosome level----------------------------

float get_chrom_disAenergy(int thread_num,int chr_homology,int chr_neighbor_homology,MatrixXf &new_chromosome,vector<MatrixXf> &tads_center_chr,vector<vector<float>> &homology_dis_matrix){

    int tad_len=chr_tad_len[chr_homology/2];
    int tad_neighbor_chr_len=chr_tad_len[chr_neighbor_homology/2];
    int indexi,indexj,indexi_haploid,indexj_haploid;
    float x,y,z,dis,dis_pre,energy_change=0;

    int homology_start = homology_tad_start[chr_homology];
    int neighbor_homology_start = homology_tad_start[chr_neighbor_homology];   

    int haploid_start = haploid_tad_start[chr_homology/2];
    int neighbor_haploid_start = haploid_tad_start[chr_neighbor_homology/2]; 

    //change neighbor_chr_tad_center 2 tad_center_vec
    vector<float> tad_center_vec(tads_center_chr[chr_neighbor_homology].data(),tads_center_chr[chr_neighbor_homology].data()+tads_center_chr[chr_neighbor_homology].size());
    vector<float> new_chromosome_vec(new_chromosome.data(),new_chromosome.data()+new_chromosome.size());

    for(int tad_index=0;tad_index<tad_len;tad_index++){
        for(int tad_neighbor_index=0;tad_neighbor_index<tad_neighbor_chr_len;tad_neighbor_index++){

            indexi=homology_start+tad_index;            
            indexj=neighbor_homology_start+tad_neighbor_index;  

            //update distance
            x=tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index]-new_chromosome_vec[tad_len*0+tad_index];   
            y=tad_center_vec[tad_neighbor_chr_len*1+tad_neighbor_index]-new_chromosome_vec[tad_len*1+tad_index];
            z=tad_center_vec[tad_neighbor_chr_len*2+tad_neighbor_index]-new_chromosome_vec[tad_len*2+tad_index];                 
            // // cout<<tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index]<<'\t'<<new_chromosome_vec[tad_len*0+indexi]<<'\n';
            dis=sqrt(x*x+y*y+z*z);
            // if(dis==0){
            //     cout<<"========================"<<'\n';
            //     cout<<tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index-1]<<'\t'<<new_chromosome_vec[tad_len*0+tad_index-1]<<'\t'<<tads_center_chr[chr_homology](tad_index-1,0)<<'\n';
            //     cout<<tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index]<<'\t'<<new_chromosome_vec[tad_len*0+tad_index]<<'\t'<<tads_center_chr[chr_homology](tad_index,0)<<'\n';
            //     cout<<tad_center_vec[tad_neighbor_chr_len*0+tad_neighbor_index+1]<<'\t'<<new_chromosome_vec[tad_len*0+tad_index+1]<<'\t'<<tads_center_chr[chr_homology](tad_index+1,0)<<'\n';
            //     cout<<thread_num<<'\t'<<chr_homology<<'\t'<<chr_neighbor_homology<<'\t'<<tad_index<<'\t'<<tad_neighbor_index<<'\n';
            
            //     cout<<tad_center_vec[tad_neighbor_chr_len*1+tad_neighbor_index]<<'\t'<<new_chromosome_vec[tad_len*1+tad_index]<<'\t'<<tads_center_chr[chr_homology](tad_index,1)<<'\n';  
            //     cout<<tad_center_vec[tad_neighbor_chr_len*2+tad_neighbor_index]<<'\t'<<new_chromosome_vec[tad_len*2+tad_index]<<'\t'<<tads_center_chr[chr_homology](tad_index,2)<<'\n';                          
            
            // }
            
            dis_pre=homology_dis_matrix[min(indexi,indexj)][max(indexi,indexj)];   
            
            homology_dis_matrix_change[thread_num][tad_index][indexj]=dis; 
            // if(indexi>tad_len){cout<<dis<<'\n';}

            // if(homology_dis_matrix_change[thread_num][tad_index][indexj]==0){
            //     cout<<1111111111111111111<<'\t'<<indexi<<'\t'<<indexj<<'\n';
            // }          

            //update energy   

            indexi_haploid=haploid_start+tad_index;
            indexj_haploid=neighbor_haploid_start+tad_neighbor_index;  
            float energy_pre=energy_change;
            // if(chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis_pre,float(1))/max(dis_pre,float(1))!=0){cout<<"error";}
            energy_change -= chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis_pre,float(1))/max(dis_pre,float(1));  
            energy_change += chr_E[min(indexi_haploid,indexj_haploid)][max(indexi_haploid,indexj_haploid)]/max(dis,float(1))/max(dis,float(1));                    

        }        
    }  

    return energy_change;              

}

float get_chrom_delta_E(int thread_num,int chr_homology,MatrixXf &new_chromosome,vector<MatrixXf> &tads_center_chr,vector<vector<float>> &homology_dis_matrix){

    int chr_another_homology,another_homology_index,tad_neighbor_chr_len,indexi,indexj,dis;
    int tad_len=chr_tad_len[chr_homology/2];
    int tad_start=homology_tad_start[chr_homology];
    float energy_change=0;

    for(int chr_neighbor=0;chr_neighbor<chr_homology_num;chr_neighbor++){ 
        if(chr_neighbor==chr_homology){continue;}
        // if(chr_neighbor/2==chr_homology/2){continue;}        
        energy_change+=get_chrom_disAenergy(thread_num,chr_homology,chr_neighbor,new_chromosome,tads_center_chr,homology_dis_matrix);                   
    }
        
    return energy_change;
}

void update_chrom_dis_matrix(int thread_num,int chr_homology){

    //update neighbor
    vector<vector<float>> change_dis;
    int indexi,indexj,index=0;
    int start_index=homology_tad_start[chr_homology];
    int end_index=homology_tad_start[chr_homology]+chr_tad_len[chr_homology/2];
    int chr_len=chr_tad_len[chr_homology/2];

    // if(chr_homology%2==1)
    // {
    //     index=chr_tad_len[chr_homology/2];
    //     start_index=homology_tad_start[chr_homology-1];
    // }
    // else{
    //     end_index=homology_tad_start[chr_homology+1]+chr_tad_len[chr_homology/2];
    // }
    
    for(int i=0;i<chr_len;i++){

        for(int j=0;j<tads_num;j++){
            indexi=start_index+i;
            indexj=j;
            if(j>=start_index && j<=end_index){
                continue;
            }            
                
            if(indexi>indexj){
                swap(indexi,indexj);
            }
            
            // homology_dis_matrix[thread_num][indexi][indexj]=homology_dis_matrix_change[thread_num][i+index][j];
            homology_dis_matrix[thread_num][indexi][indexj]=homology_dis_matrix_change[thread_num][i][j];
            // if(homology_dis_matrix_change[thread_num][i+index][j]==0){
            //     cout<<start_index<<'\t'<<end_index<<'\t'<<i<<'\t'<<(i+start_index)<<'\t'<<indexi<<'\t'<<indexj<<'\n';
            // }

        }
    }  
}

void update_chrom_neighbor_contact(int homology,int thread_num,int m){
    int m_pre=homology_m[thread_num][homology];
    int chr_len=chr_tad_len[homology/2];
    int homology_start=homology_tad_start[homology];
    int homology_end=homology_start+chr_len;
    float dis;
    float contact;

    int chr_homology_index=homology%2;
    // if(chr_homology_index==0){homology_end=homology_tad_start[homology+1]+chr_len;}
    // else{homology_start=homology_tad_start[homology-1];}

    for(int i=homology_start;i<homology_end;i++){
        for(int j=0;j<tads_num;j++){
            // if(j>=homology_start && j<homology_end){continue;}
            dis=min(homology_dis_matrix[thread_num][min(i,j)][max(i,j)],float(40));
            contact=float(m-m_pre)/pow(max(dis,float(1)),2.3);
            homology_contact_matrix[thread_num][min(i,j)][max(i,j)]+=contact;
            // if(homology_contact_matrix[thread_num][min(i,j)][max(i,j)]<0){cout<<i<<'\t'<<j<<'\n';}
        }         
    }   
    // cout<<contact<<'\t'<<homology_start<<'\t'<<homology_contact_matrix[thread_num][max(homology_start-10,0)][homology_start+10]<<'\n';       

    homology_m[thread_num][homology]=m;
}


void normalize_chrom() {   

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
            // cout<<sum<<'\t';
            // if(haploid_contact_matrix[i][j]<0){cout<<i<<'\t'<<j<<'\t'<<haploid_contact_matrix[i][j]<<'\n';}
        }
    }

    for(int i=0;i<chr_homology_num/2;i++){
        int chr_start=haploid_tad_start[i];
        int chr_end=chr_tad_len[i]+chr_start;

        for(int j=chr_start;j<chr_end;j++){
            for(int k=j+1;k<chr_end;k++){
                sum_intra+=haploid_contact_matrix[j][k];
                sum_intra_refer+=reference_genome_tads_contact[j][k];

                haploid_contact_matrix[j][k]=0;
                haploid_contact_matrix[k][j]=0;

                // reference_genome_tads_contact[j][k]=0;
                // reference_genome_tads_contact[k][j]=0;                
            }
        }
    }

    sum_inter_refer=sum_refer-sum_intra_refer;
    sum_inter=sum-sum_intra;
    float multiple_inter=sum_inter_refer/sum_inter;

    cout<<"sum:"<<sum_inter_refer<<'\t'<<sum_refer<<'\t'<<sum<<'\n';

    //normalize
    for(int i=0;i<length;i++){
        for(int j=i+1;j<length;j++){
            haploid_contact_matrix[i][j]=haploid_contact_matrix[i][j]*multiple_inter;
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

void update_chrom_energies1() {
    //update energy & combine real contact and refer contact
    float checksum = 0;
    float loss=0;
    int length=tads_num/2;
    float margin;

    for (int i = 0; i < length; i++) {

        //update energy
        for (int j = i+1; j < length; j++) {
            margin=haploid_contact_matrix[i][j]-reference_genome_tads_contact[i][j];         
            // chr_E[i][j] += learning_rate * sqrt(1/(max(reference_genome_tads_contact[i][j],float(0.0001))))*margin;
            // chr_E[j][i] = chr_E[i][j];
            loss+=fabs(margin);        

            //combine real contact and refer contact
            haploid_contact_matrix[j][i]=reference_genome_tads_contact[i][j];
        }
    }

    cout<<"loss"<<"\t"<<loss<<"\n";
}

void update_chrom_energies2() {
   //update energy & combine real contact and refer contact
    float checksum = 0;
    float loss=0;
    int length=tads_num/2;
    float margin;

    for (int i = 0; i < length; i++) {

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

    //mask
    for(int i=0;i<chr_homology_num/2;i++){
        int chr_start=haploid_tad_start[i];
        int chr_end=chr_tad_len[i]+chr_start;

        for(int j=chr_start;j<chr_end;j++){
            for(int k=j+1;k<chr_end;k++){

                chr_E[j][k]=0;
                chr_E[k][j]=0;         
            }
        }
    }    


}

void update_chrom_energies() {
    //update energy & combine real contact and refer contact
    float checksum = 0;
    float loss=0,refer_contact=0;
    int length=tads_num/2;
    float margin;

    //combine real contact and refer contact
    for (int i = 0; i < length; i++) {
        for (int j = i+1; j < length; j++) {
            haploid_contact_matrix[j][i]=reference_genome_tads_contact[i][j];
        }
    }

    //mask
    for(int i=0;i<chr_homology_num/2;i++){
        int chr_start=haploid_tad_start[i];
        int chr_end=chr_tad_len[i]+chr_start;

        for(int j=chr_start;j<chr_end;j++){
            for(int k=j+1;k<chr_end;k++){

                haploid_contact_matrix[j][k]=0;
                haploid_contact_matrix[k][j]=0;         
            }
        }
    }
            
    //update energy
    for (int i = 0; i < length; i++) {
        for (int j = i+1; j < length; j++) {

            refer_contact=haploid_contact_matrix[j][i];
            margin=haploid_contact_matrix[i][j]-refer_contact;         
            
            chr_E[i][j] += learning_rate * sqrt(1/(max(refer_contact,float(0.00001))))*margin;
            chr_E[j][i] = chr_E[i][j];
            loss+=fabs(margin);        
        }
    }    

    cout<<"loss"<<"\t"<<loss<<"\n";
     
}




#endif //CHROM_ENERGY_CHANGES_H