#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <boost/functional/hash.hpp>

 
using namespace Eigen;
using namespace std;

struct pair_hash {
    size_t operator()(const pair<int, int>& pair) const noexcept {
        size_t seed = 0;
        boost::hash_combine(seed, pair.first);
        boost::hash_combine(seed, pair.second);//把pair.first 和pair.second combine映射成一个哈希值
        return seed;
    }
};

int main()
{
    unordered_map<pair<int, int>, double, pair_hash> dict;

    for(int i=0;i<1000;i++){
        for(int j=0;j<1000;j++){
            dict[{i,j}]=i+j;
        }
    }

    typedef Eigen::Triplet<double> T;
    std::vector<T> c;
    auto start = chrono::high_resolution_clock::now();
    int a;
    for(auto &tmp:dict){
        dict[tmp.first]+=1;
        cout<<dict[tmp.first]<<'\t';
        c.push_back(T(tmp.first.first,tmp.first.second,dict[tmp.first]));
    }
  
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<float> elapsed1=end-start;
    cout<<elapsed1.count()<<'\n';    
    Eigen::SparseMatrix<double> mat(10000,10000);
    mat.setFromTriplets(c.begin(),c.end());
    // cout<<mat;

    // SparseMatrix<float> mat(10,10),spMat;
    // // A = MatrixXf::Random(5,5).sparseView(); 

    // // int k=2;
    // // vector<int> k_list={1,2,3};

    // // for(auto i=0;i<A.cols();i++){
    // //     A.prune([&k_list](int i, int j, float) { return i!=100 & j!=100;  }); 
    // // }
    // A.setZero();
    // // A.coeffRef(0,1)=5.0;
    // int i=0;int j=1;
    // // A.insert(i,j)=3;
    // auto start = chrono::high_resolution_clock::now();
    // for(auto i=0;i<A.cols();i++){
    //     for(auto j=0;j<A.rows();j++){
    //         A.coeffRef(j,i);        
    //     }
    // }    
    // auto end = chrono::high_resolution_clock::now();
    // chrono::duration<float> elapsed1=end-start;
    // cout<<elapsed1.count()<<'\n';

    // cout<<"finish"<<'\t';

    return 0;
}