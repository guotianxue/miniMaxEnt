#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

 
using namespace Eigen;
using namespace std;
int main(int argc, char** argv)
{
    SparseMatrix<float> mat(100,100),spMat;

    mat.resize(10,10); 
    SparseMatrix<float>  sm1;
    sm1=mat;
    sm1.coeffRef(2,1) = 0.5;
    MatrixXf dmat;
    dmat = MatrixXf(sm1);
    spMat = dmat.sparseView();
    mat.setZero(); 
    cout<<mat;
    // ofstream Sout("test.txt");
    // Sout << dmat;
    // Sout.close();      


    return 0;
}