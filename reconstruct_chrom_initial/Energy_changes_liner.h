#ifndef INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H
#define INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
//#include "Functions.h"
#include "global.h"
using namespace std;
//----------------calculate logic------------
bool is_equal_loc(vector<int> loc1,vector<int> loc2){
    int len1=loc1.size();
    int len2=loc2.size();
    if(len1!=len2){return false;}
    for(int i=0;i<len1;i++){
        if(loc1[i]!=loc2[i]){
            return false;
        }
    }
    return true;
}

vector<int> add_loc(vector<int> loc1,vector<int> loc2,int k1,int k2){
    int len1=loc1.size();
    int len2=loc2.size();
    vector<int> res(len1);

    for(int i=0;i<len1;i++){
        res[i]=k1*loc1[i]+k2*loc2[i];
    }
    return res;
}


//-----------------------chrom level---------------------------
float delta_E_chrom(vector<vector<int>> &polymer,vector<int> pre_move,vector<int> prop_move,int chromosome_index,int homology){
    //update neighbor

    float dis;float dis_pre;vector<int> margin;vector<int> margin_pre;
    float energy_change=0;
    int r_bound;
    int index=chromosome_index*2+homology;
    for(int i=0;i<polymer.size();++i){
        if(i!=chromosome_index*2 && i!=chromosome_index*2+1){
            margin_pre=add_loc(polymer[i],polymer[index],1,-1);
            margin=add_loc(polymer[i],prop_move,1,-1);
            dis=sqrt(pow(margin[0],2)+pow(margin[1],2)+pow(margin[2],2)); 
            dis_pre=sqrt(pow(margin_pre[0],2)+pow(margin_pre[1],2)+pow(margin_pre[2],2)); 

            r_bound=chr_r_list[i/2]+chr_r_list[chromosome_index/2];
            
            if(beta>500){
                energy_change+=min(dis_pre-r_bound,float(0))*(beta-500)/2000;
                energy_change-=min(dis-r_bound,float(0))*(beta-500)/2000;                
            }

            energy_change -= Interaction_E_chrom(chromosome_index,i/2)*(pow(dis_pre,-1));//*(pow(dis_pre,-1))  ;
            energy_change += Interaction_E_chrom(chromosome_index,i/2)*(pow(dis,-1));//*(pow(dis,-1));
        }
    }

    return energy_change;

}

void update_chrom_contact(vector<vector<int>> &polymer, int thread_num){

    //Throw away old contacts, at the same time update contact frequency map 更新旧的site位点接触的contact
    int chromosome_index;
    vector<int> margin;
    float dis;

    for(int i=0;i<polymer.size();++i){
        chromosome_index=i/2;
        for(int j=i+1;j<polymer.size();j++){
            if(j/2!=chromosome_index ){
                margin=add_loc(polymer[i],polymer[j],1,-1);
                dis=sqrt(pow(margin[0],2)+pow(margin[1],2)+pow(margin[2],2));
                if(dis==0){
                    chrom_total_contacts[thread_num](min(chromosome_index,j/2),max(chromosome_index,j/2)) += 1; 
                }
                else{
                    chrom_total_contacts[thread_num](min(chromosome_index,j/2),max(chromosome_index,j/2)) += pow(dis,-1.5); 
                }
            }
        }
    }

}


void update_energy_chrom(Eigen::MatrixXd &refer) {

    float checksum = 0;
    float loss=0;
    Eigen::MatrixXd margin;

    margin=chrom_final_contacts- refer;
    for (int i = 0; i < chr_num; i++) {
        for (int j = i+1; j < chr_num; j++) {

            Interaction_E_chrom(i,j) = learning_rate * sqrt(1/(max(refer(i,j),double(0.0001))))*(chrom_final_contacts(i,j)- refer(i,j));
            Interaction_E_chrom(j,i)=Interaction_E_chrom(i,j);
            loss+=fabs(margin(i,j));   
            // cout<< learning_rate * sqrt(1/(max(refer(i,j),double(0.0001))))*(chrom_final_contacts(i,j)- refer(i,j))<<'\t';                   
        }
    }
    cout<<"loss"<<"\t"<<loss<<"\n";

    // //calculates the shift of all energies, imposed to ensure a MaxEnt solution is found for the contact frequency scale
    // float shift=2*((Interaction_E_chrom.cwiseProduct(refer)) / chr_num).sum();
    // cout << "Shift: " << shift << endl;

    // //  //update E
    // Eigen::MatrixXd matrix_multiple=Eigen::MatrixXd::Ones(chr_num,chr_num);
    // Interaction_E_chrom-=shift*matrix_multiple;

}

//Normalizes model contact frequencies to allow for comparison with experimental data
void normalize() {

    float sum = chrom_final_contacts.sum();
    float multiply= chr_num/(2*sum);

    chrom_final_contacts*= multiply;
    cout<<"nor_check"<<"\t"<<chrom_final_contacts.sum()<<"\n";
}

void normalize_reference(Eigen::MatrixXd &refer) {

    float sum=refer.sum();
    refer*=chr_num/(2*sum);

    cout<<"refer_check"<<"\t"<<refer.sum()<<"\n";
}

#endif //INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H