//
// Created by Joris on 13/07/2018.
//
#ifndef INVERSE_3D_NEW_GOBAL_H
#define INVERSE_3D_NEW_GOBAL_H


#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include <string>

#include <boost/functional/hash.hpp>
using namespace std;
using namespace Eigen;

extern string chr;
extern int chr_len;
extern int pol_length;
extern int tad_bin_num;
extern const int refer_length;
extern const int mid_pol_length;
extern long int mc_moves_start;
extern long int mc_moves;
extern long int burn_in_time;
extern int number_of_threads;
extern const int cap;
extern float learning_rate;

extern int mask_bins;
extern vector<int> mask_bin_list;

//neigbor parameter
extern int threshold;//if bin dis >=200 ,deal with neighbor in 3d
extern int neighbor_dis;//if dis>neighbor_dis,regard as a neighbor contact



// extern static mt19937_64 gen;
extern uniform_real_distribution<float> unif;
extern uniform_int_distribution<int> unimove;
extern uniform_int_distribution<int> unisite;
extern uniform_int_distribution<int> unipath;//in,out,(x)left,right(y),up,down(z)
extern uniform_int_distribution<int> leftOright;
extern uniform_int_distribution<int> uniradius;

extern const int length_cylinder;
extern const int cap_length;
extern const int midway;

extern vector<vector<vector<int>>> polymer;
extern vector<vector<vector<int>>> polymer_bin;

extern vector<pair<int, int>> tads;

//record polymer site index (vector) ------> (site index) : (polymer position site)
extern vector<vector<int>> site_index2position;

//record polymer site index (vector) ------> (polymer position site) : (site index)
extern vector<vector<int>> position2site_index;

//record site available flag     0 means no site   1 means have site
extern vector<vector<int>> site_available;


extern vector<vector<vector<float>>> total_contacts;
extern vector<vector<int>> contacts_list;
extern vector<vector<int>> prop_contacts_list;

//record the gene site position 
extern vector<int> site_position ;
//record the number of points behind each genome site
extern unordered_map<int, int> site_behind_num;


extern bool boundary_cond;
extern bool orient;
extern const float diameter;
extern int radius;
extern const int length_cell;
extern const int midway;

struct pair_hash {
    size_t operator()(const pair<int, int>& pair) const noexcept {
        size_t seed = 0;
        boost::hash_combine(seed, pair.first);
        boost::hash_combine(seed, pair.second);//把pair.first 和pair.second combine映射成一个哈希值
        return seed;
    }
};
struct vec_hash {
    size_t operator()(const vector<int> vec) const noexcept {
        return boost::hash_range(vec.cbegin(), vec.cend());//把数组映射成一个哈希值
    }
};

struct NeighborDefaulted {
    float dis=0;
    int m=0;
};

//all contact initial
extern unordered_map<pair<int, int>, float, pair_hash> reference_contacts_num;
extern unordered_map<int, vector<int>> reference_contacts_pair;
extern SparseMatrix<float> reference_contact_s;

//every tad initial
extern vector<unordered_map<int, vector<int>>> reference_tad_contacts_pair;//(tad_numbers,pair,connect_pair_vector)
extern vector<unordered_map<pair<int, int>, float, pair_hash>> reference_tad_contacts_num;//(tad_numbers,pair,numOfContacts)

//initial mask
// extern Eigen::MatrixXf mask;
extern SparseMatrix<float> mask_s;

//one tad initial
extern SparseMatrix<float> reference_one_tad_contact_s;    
extern Eigen::MatrixXf refer_msk_contact;  

extern vector<unordered_map<pair<int, int>, float, pair_hash>> one_tad_contact;
extern vector<unordered_map<vector<int>, vector<int>, vec_hash>> one_tad_locations;

extern vector<SparseMatrix<float>> tad_total_contacts_s;
extern SparseMatrix<float> tad_final_contacts_s;
extern vector<SparseMatrix<float>> tad_neighbor_contacts_s;

//initial neighbor
extern vector<unordered_map<pair<int, int>, float, pair_hash>> tad_total_contacts_neighbor;
extern vector<unordered_map<pair<int, int>, float, pair_hash>> neighbor_contact_dis_hash;
extern vector<unordered_map<pair<int, int>, int, pair_hash>> neighbor_contact_m_hash;
extern vector<unordered_map<pair<int, int>, NeighborDefaulted, pair_hash>> neighbor_contact_data;

extern vector<unordered_map<int, vector<int>>> tad_neighbor_contact;
extern vector<unordered_map<int, vector<int>>> tad_neighbor_locations;
extern vector<vector<int>> tad_neighbor_contact_new;

extern vector<unordered_map<pair<int, int>, float, pair_hash>> neighbor_contact_dis;

//energy initial 
extern SparseMatrix<float> Interaction_E_tad_s;//tad中相互作用的能量
extern Eigen::MatrixXf Interaction_E_tad;
extern unordered_map<pair<int, int>, float, pair_hash> Interaction_E_tad_1;

//initial  multiple matrix
extern Eigen::MatrixXf matrix_multiple;
extern Eigen::MatrixXf deal_multiple;

extern Eigen::MatrixXf  ref_final;

//test
extern int reduce_count,add_count,kink_count,loop_count,skip_count;

extern unordered_map<pair<int, int>, float, pair_hash> refer_hash;

#endif //INVERSE_3D_NEW_GOBAL_H

