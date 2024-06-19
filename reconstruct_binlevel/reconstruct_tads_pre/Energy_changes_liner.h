#ifndef INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H
#define INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
//#include "Functions.h"
#include "global.h"
using namespace std;


//--------------------calculate delta E------------------------
float delta_E_other(vector<vector<int>> &polymer ,vector<int> &site_index2position, vector<int> &site_available, int site, int pol_length, vector<int> prop_move1,int thread_num){
    //求kink和loop move的detaE  改变site位置
    if (site%6!=0){
        // not in bin site
        return 0;
    }

    float energy_change = 0;
    float dis;
    //弄一个energy有值的hashmap啦

    if (one_tad_locations[thread_num].find(prop_move1)!=one_tad_locations[thread_num].end()){//如果prop_move1 (新的第二个位点) 在location中
        for (auto elem : one_tad_locations[thread_num].find(prop_move1)->second ) {
            if (elem%6==0 ){
                energy_change += Interaction_E_tad(min(elem/6,site/6),max(elem/6,site/6));
            }
        }
    }

    for (auto elem : one_tad_locations[thread_num].find(polymer[site_index2position[site]])->second ) {
        if (elem%6==0 and elem != site){      
            energy_change -= Interaction_E_tad(min(elem/6,site/6),max(elem/6,site/6));//因为原来的第二个位点(site+1)被改变了，应减去location原来每个点与点site+1接触的能量
        }
    }

    for(auto elem : tad_neighbor_contact_new[thread_num]){
        dis=neighbor_contact_data[thread_num][{min(elem,site),max(elem,site)}].dis;
        if(dis!=0){
            energy_change +=  Interaction_E_tad(min(elem/6,site/6),max(elem/6,site/6))*pow(dis,-1)*0.5;//*pow(max(-dis,float(1)),1) ;// /pol_length*(max(site,elem[0])-min(site,elem[0]));
        }

    }

    if(tad_neighbor_contact[thread_num].find(site)!=tad_neighbor_contact[thread_num].end()){
        for(auto elem : tad_neighbor_contact[thread_num][site]){
            dis=neighbor_contact_data[thread_num][{min(elem,site),max(elem,site)}].dis;   
            if(dis!=0) {
                energy_change-= Interaction_E_tad(min(elem/6,site/6),max(elem/6,site/6))*pow(dis,-1)*0.5;//*pow(max(-dis,float(1)),1) ;// /pol_length*(max(site,elem[0])-min(site,elem[0]));
            }    
        }
    }

    return energy_change;
}

//--------------------------update contact-----------------------------
void update_contact_location(vector<vector<int>> &polymer, int site, int pol_length, vector<int> prop_move1,vector<int> pre_site,int thread_num,int m){

    //Throw away old contacts, at the same time update contact frequency map 更新旧的site位点接触的contact

    if ((site%6)!=0) return;

    for (auto elem : one_tad_locations[thread_num][pre_site]){
        if ( abs(elem-site)>mask_bins*6){
            tad_total_contacts_s[thread_num].coeffRef(min(elem,site)/6,max(elem,(site))/6) += m - one_tad_contact[thread_num][{min(elem,(site)),max(elem,site)}];
            // cout<<tad_total_contacts_s[thread_num].coeffRef((min(elem,(site)))/6,(max(elem,(site)))/6)<<'\n';
            one_tad_contact[thread_num].erase({min(elem,(site)),max(elem,(site))});
     
        }
    }
    //put in new contacts
    if (one_tad_locations[thread_num].find(prop_move1)!=one_tad_locations[thread_num].end()){//新site与location中位置重合
        for (auto elem : one_tad_locations[thread_num].find(prop_move1)->second ) {
            if ( abs(elem-site)>mask_bins*6){
                one_tad_contact[thread_num][{min(elem,site),max(elem,site)}] = m;
                // cout<<"k_add"<<"\t"<<elem<<"\t"<<site<<"\t"<<one_tad_contact[thread_num][{min(elem,site),max(elem,site)}]<<"\n";            
            }
        }
    }

    //update hash map locations
    vector<int> first_monomer = pre_site;
    if (one_tad_locations[thread_num][first_monomer].size() ==1){//first_monomer位置只有site一个位点，更新之后直接删
        one_tad_locations[thread_num].erase(first_monomer);
    }
    else {
        one_tad_locations[thread_num][first_monomer].erase(find(one_tad_locations[thread_num][first_monomer].begin(), one_tad_locations[thread_num][first_monomer].end(),site));
        //擦除first_monomer的key对应vector中的site位点
    }
    //将location更新，改为新的site的坐标
    if (one_tad_locations[thread_num].find(prop_move1) != one_tad_locations[thread_num].end()){
        one_tad_locations[thread_num][prop_move1].emplace_back(site);
    }
    else {
        one_tad_locations[thread_num][prop_move1] = {site};
    }
}

void update_contact_location_neighbor( int site, int thread_num,int m){

    //Throw away old contacts, at the same time update contact frequency map 
    if (site%6!=0) {return;}

    float dis;float dis_new;int mc_move;int index;
    for (auto elem:tad_neighbor_contact_new[thread_num]){

        dis=neighbor_contact_data[thread_num][{min(site,elem),max(site,elem)}].dis;        
        if(dis==0 && abs(site-elem)>mask_bins*6){
            dis_new=neighbor_contact_dis[thread_num][{min(site,elem),max(site,elem)}];

            tad_neighbor_contact[thread_num][elem].emplace_back(site);
            //update
            neighbor_contact_data[thread_num][{min(site,elem),max(site,elem)}].m=m;
            neighbor_contact_data[thread_num][{min(site,elem),max(site,elem)}].dis=dis_new;            
        }      
    }        

    for(auto elem : tad_neighbor_contact[thread_num][site]){
        if(abs(site-elem)>mask_bins*6){

            dis=neighbor_contact_data[thread_num][{min(site,elem),max(site,elem)}].dis;
            dis_new=neighbor_contact_dis[thread_num][{min(site,elem),max(site,elem)}];
            mc_move=neighbor_contact_data[thread_num][{min(site,elem),max(site,elem)}].m;

            tad_total_contacts_neighbor[thread_num][{(min(elem,site))/6,(max(elem,site))/6}]+= float(m - mc_move)*pow(max(dis,float(1)),-3);
            
            //update
            if(dis_new==0 ){
                vector<int> &elem_contact=tad_neighbor_contact[thread_num][elem];
                elem_contact.erase(find(elem_contact.begin(), elem_contact.end(),site));

            }

            neighbor_contact_data[thread_num][{min(site,elem),max(site,elem)}].m=m;
            neighbor_contact_data[thread_num][{min(site,elem),max(site,elem)}].dis=dis_new;            
        }
    }    
    
    tad_neighbor_contact[thread_num].erase(site);

    //put in new contacts,update contact and locations
    if (tad_neighbor_contact_new[thread_num].empty()==0){
        tad_neighbor_contact[thread_num][site]=tad_neighbor_contact_new[thread_num];
    } 
}


//--------------------------------------------------------------------------------------------------------
//after each forward run, the polymer interaction energies are updated according to the pairwise difference
// between model contact frequencies and experimental contact frequencies

void update_energies(SparseMatrix<float> &refer) {

    float checksum = 0;
    float loss=0;
    SparseMatrix<float> margin;

    cout<<"update"<<'\n';

    margin=tad_final_contacts_s- refer;

    cout<<111<<'\n';

    int i,j;
    float refer_data,margin_data;
    for (int k=0; k<margin.outerSize(); ++k){
        for (SparseMatrix<float>::InnerIterator it(margin,k); it; ++it)
        {
            loss+=fabs(it.value());
            i=it.row();
            j=it.col();            

            refer_data=refer_hash[{it.row(),it.col()}];
            margin_data=it.value();
            Interaction_E_tad(i,j) += learning_rate * sqrt(1/(max(refer_data,float(0.0001))))*margin_data;   
            // Interaction_E_tad(i,j) += learning_rate * sqrt(1/(max(refer_data,float(0.0001))))*margin_data;   
            // Interaction_E_tad_1[{i,j}]+=learning_rate * sqrt(1/(max(refer_data,float(0.0001))))*margin_data;
        }
    }    
     
    cout<<"loss"<<"\t"<<loss<<"\n";

    //update energy hashmap

    // //calculates the shift of all energies, imposed to ensure a MaxEnt solution is found for the contact frequency scale
    // float shift=(2*(Interaction_E_tad_s.cwiseProduct(refer)) / tad_bin_num).sum();
    // cout << "Shift: " << shift << endl;

    //  //update E
    // Interaction_E_tad-=shift*matrix_multiple;
    
    cout<<"energy_sum:"<<Interaction_E_tad.sum()<<'\n';   
}


//Normalizes model contact frequencies to allow for comparison with experimental data


#endif //INVERSE_3D_NEW_ENERGY_CHANGES_LINER_H