#ifndef INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H
#define INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
//#include "Functions.h"
#include "global.h"
using namespace std;
using namespace Eigen;

//============================= tads calculate ============================
//--------------------calculate delta E------------------------
float delta_E_tads( MatrixXi tad ,MatrixXi change_tad ,MatrixXf tads_center_change, int tad_start,int thread_num){
    //求kink和loop move的detaE  改变site位置

    float energy_change = 0;

    vector<int> monomer;
    vector<int> change_monomer;
    int start=tads_bound[tad_start][0];
    int end=tads_bound[tads_bound.size()-1][1];
    int loc;
    int elem_tad_index;

    //update neighbor
    Eigen::MatrixXf tads_center_new=Eigen::MatrixXf::Zero(tads_bound.size(),3);
    tads_center_new=tads_center[thread_num];
    tads_center_new.bottomRows(tads_bound.size()-tad_start)=tads_center_change;

    float dis;float dis_pre;Eigen::MatrixXf margin;
    for(int i=tad_start;i<tads_bound.size();++i){
        for(int j=0;j<tads_bound.size();++j){
            if(i!=j){//if(abs(i-j)<=200 && i!=j)
                margin=tads_center_new.row(i)-tads_center_new.row(j);           
                dis=sqrt((margin.cwiseProduct(margin)).sum()); 
                dis_pre=tads_center_dis[thread_num](i,j);
                if (dis_pre<=10){
                    energy_change -= Interaction_E_tads(i,j)*dis_pre;//*(pow(dis_pre,-0.33))  ;
                }
                if (dis<=10){
                    energy_change += Interaction_E_tads(i,j)*dis;//*(pow(dis,-0.33));
                }
                // cout<<dis_pre<<'\t'<<dis<<'\n';
            }
        }
    }    

    return energy_change;
}


//------------------------update tad contact-----------------------------
void update_tads_location(MatrixXi tad ,MatrixXi change_tad ,int tad_start_index,int thread_num,int m){

    vector<int> monomer;
    vector<int> change_monomer;
    int start;int end;int loc;float dis;
    int first_index=tads_bound[0][0]; 
    int first=tads_bound[tad_start_index][0];   
    int last=tads_bound[tads_bound.size()-1][1];   
    //Throw away old contacts, at the same time update contact frequency map 更新旧的site位点接触的contact
 
    //update neighbor
    for(int i=0;i<tads_bound.size();i++){
        for(int j=0;j<tads_bound.size();j++){
            dis=tads_center_dis[thread_num](i,j);
            if(i!=j && dis<=10){
                tads_total_contacts[thread_num](min(i,j),max(i,j)) +=min(pow(dis,-0.33),10.00);        
            }
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

    //update neighbor
    for(int i=0;i<tads_bound.size();i++){
        for(int j=0;j<tads_bound.size();j++){
            dis=tads_center_dis[thread_num](i,j);
            if(i!=j && dis<=10){
                tads_total_contacts[thread_num](min(i,j),max(i,j)) +=min(pow(dis,-0.33),10.00);        
            }
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

    int length=tads_bound.size();
    float sum = tads_final_contacts.sum();
    float check1=0;
    float multiply= float(length)/(2*sum);

    tads_final_contacts*= multiply;// pol_length/4     /(8*sum)=/(4*4*sum)*2
    check1+=tads_final_contacts.sum();         


    cout<<"nor_check"<<"\t"<<check1<<"\n";
}

void normalize_tads_reference() {

    // ofstream no_refer_contact;
    int length=tads_bound.size();
    float sum = reference_tads_contact.sum();
    reference_tads_contact*=length/(2*sum);
    
    cout<<"refer_check"<<"\t"<<reference_tads_contact.sum()<<"\n";
}

void update_tads_energies() {
    float checksum = 0;
    float loss=0;
    int length=tads_bound.size();
    Eigen::MatrixXf margin;
    Eigen::MatrixXf One=Eigen::MatrixXf::Ones(length,length)-Eigen::MatrixXf::Identity(length,length) ;

    margin=tads_final_contacts- reference_tads_contact;
    for (int i = 0; i < length; i++) {
        for (int j = i+1; j < length; j++) {
            Interaction_E_tads(i,j) += learning_rate * sqrt(1/(max(reference_tads_contact(i,j),float(0.001))))*(tads_final_contacts(i,j)- reference_tads_contact(i,j));
            Interaction_E_tads(j,i) = Interaction_E_tads(i,j);
            loss+=fabs(margin(i,j));
            checksum += (tads_final_contacts(i,j))*length;//checksum 512 (32*32/2)                        
        }
    }
    cout<<"loss"<<"\t"<<loss<<"\n";

    //calculates the shift of all energies, imposed to ensure a MaxEnt solution is found for the contact frequency scale
    
    float shift=(2*(Interaction_E_tads.cwiseProduct(reference_tads_contact)) / length).sum();

    cout << "Shift: " << shift << endl;
    //update E
    Interaction_E_tads-=shift*One;
    cout<<"lala:"<<(Interaction_E_tads).cwiseProduct(reference_tads_contact+reference_tads_contact.transpose()).sum()<<'\n';
    cout<<"energy_sum:"<<Interaction_E_tads.sum()<<'\n';   
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


#endif //INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H