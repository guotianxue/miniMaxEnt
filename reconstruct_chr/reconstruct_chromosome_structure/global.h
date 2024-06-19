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
extern int chr_homology_num;
extern int pol_length;
extern int tad_bin_num;
extern int tads_num;
extern int mark;
extern int reconstruct_chrom_mark;
extern const int refer_length;
extern const int mid_pol_length;
extern long int mc_moves_start;
extern long int mc_moves;
extern long int burn_in_time;
extern int number_of_threads;
extern const int cap;
extern float learning_rate;

extern vector<int> is_move_active;

//neigbor parameter
extern int threshold;//if bin dis >=200 ,deal with neighbor in 3d
extern int neighbor_dis;//if dis>neighbor_dis,regard as a neighbor contact
extern float neighbor_energy_weight; 
extern float neighbor_contact_weight;


// extern static mt19937_64 gen;
extern uniform_real_distribution<float> unif;
extern uniform_int_distribution<int> unimove;
extern uniform_int_distribution<int> unisite;
extern uniform_int_distribution<int> unipath;//in,out,(x)left,right(y),up,down(z)
extern uniform_int_distribution<int> leftOright;
extern uniform_int_distribution<int> uniradius;

extern uniform_int_distribution<int> uniTranslate;

extern const int length_cylinder;
extern const int cap_length;
extern const int midway;

extern bool boundary_cond;
extern bool orient;
extern const float diameter;
extern int radius;
extern const int length_cell;
extern const int midway;
// extern vector<vector<int>> tads_bound;

extern vector<vector<int>> tads_bound_genome;
extern vector<vector<int>> tads_bound_chr;

extern vector<vector<vector<float>>> homology_dis_matrix;
extern vector<vector<vector<float>>> homology_dis_matrix_change;

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

//=============reconstruct all tads structure global variable============

extern float alpha;
extern uniform_int_distribution<int> uniThread;
extern Eigen::MatrixXf Interaction_E_tads;
extern Eigen::SparseMatrix<float> reference_tads_contact_s;
extern SparseMatrix<float> reference_genome_tads_contact_s;
extern SparseMatrix<float> reference_tads_contact_s;
extern vector<vector<float>> reference_genome_tads_contact;

extern vector<vector<int>> tads_bound;
extern Eigen::MatrixXi tad_match;

// extern unordered_map<int, int>  homology_tad_start;//chr_index:homology_start_index
// extern unordered_map<int, int>  haploid_tad_start;

extern vector<int>  haploid_tad_start;
extern vector<int>  homology_tad_start;
extern vector<int> chr_tad_len;
// extern unordered_map<int, int>  chr_tad_len;
extern vector<vector<int>> chr_tad2homology_tad;
extern vector<vector<int>> homology_m;
extern vector<int> thread_m;

extern vector<vector<vector<float>>> homology_contact_matrix;
extern vector<vector<float>> haploid_contact_matrix;

//initial contact and locations
extern vector<SparseMatrix<float>> tads_total_contacts_s;

extern SparseMatrix<float> tads_final_contacts_s;

extern vector<unordered_map<pair<int, int>, float, pair_hash>> tads_contact;
extern vector<unordered_map<vector<int>, vector<int>, vec_hash>> tads_location;
extern vector<vector<Eigen::MatrixXi>> tads;
extern vector<vector<Eigen::MatrixXi>> connected_tads;
extern vector<Eigen::MatrixXi> connected_structure;
extern vector<int> is_move_active;
extern int tad_change_flag;
extern unordered_map<int, int>  tad_bin_dict;
extern vector<Eigen::MatrixXf> tads_center;
extern vector<Eigen::MatrixXf> tads_center_m;
extern vector<Eigen::MatrixXf> tads_center_dis;
extern vector<unique_lock<mutex>*> locks;
extern std::condition_variable cv;

//==============reconstruct all bins of tads===============
extern Eigen::MatrixXd Interaction_E;
extern vector<vector<float>> chr_E;
extern Eigen::MatrixXd refer_contact;

extern Eigen::MatrixXd final_contacts;
extern Eigen::MatrixXf mask;
extern SparseMatrix<float> mask_s;

#endif //INVERSE_3D_NEW_GOBAL_H

