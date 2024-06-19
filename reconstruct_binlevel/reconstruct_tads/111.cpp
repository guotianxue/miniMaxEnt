#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>
#include <unordered_map>

 
using namespace Eigen;
using namespace std;


vector<thread> threads(10);
volatile unordered_map<int,int>  *tmp;

void run(){
    for (int i=0;i<1000;i++){
        // int a=1;
    }
}

int main()
{

    // for (auto l = 0; l < 10; l++) {
    //     threads[l] = thread(run);   
    // }
    // for (auto &&l : threads) {
    //     l.join();
    // }
    SparseMatrix<float> A(1000,1000);
    A.coeffRef(1,7)=3;
    Eigen::saveMarket(A, "./text.mtx");
    // Eigen::saveMarket(A, "filename_SPD.mtx", Eigen::Symmetric); // if A is symmetric-positive-definite
    // Eigen::saveMarketVector(B, "filename_b.mtx");


    return 0;
}