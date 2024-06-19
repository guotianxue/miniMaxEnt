#include <iostream>
#include <fstream>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <random>
#include <string>
#include <thread>
#include "reconstruct.h"
using namespace std;

std::vector<std::vector<float>> refer(100, std::vector<float>(100, 1));
std::vector<std::vector<float>> real(100, std::vector<float>(100, 1));
int main()
{ 
    float c=get_corr(refer,refer);
}