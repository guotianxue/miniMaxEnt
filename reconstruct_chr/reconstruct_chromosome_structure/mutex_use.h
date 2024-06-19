#ifndef Mutex_Initialize_h
#define Mutex_Initialize_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <random>
#include <string>
#include <thread>
#include <mutex>
using namespace std;

vector<mutex> mutexs(48);
vector<unique_lock<mutex>*> locks(48);

#endif 