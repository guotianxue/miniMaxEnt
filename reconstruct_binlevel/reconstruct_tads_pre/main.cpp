#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>  // getopt 函数所在头文件 
#include <getopt.h>  // getopt_long 函数所在头文件
#include <assert.h>
#include "reconstruct.h"

using namespace std;

static struct option long_options[] = {
    {"help", no_argument,       0,   'h'},
    {"chr",  required_argument, 0,   'H'},
    {"radius",  required_argument, 0,   'r'},
    {"learning_rate",  required_argument, 0,   'l'},
    {"number_of_threads",  required_argument, 0,   't'},
    {"update_steps",  required_argument, 0,   's'},
    {"burn_in_time",  required_argument, 0,   'b'},
    {"mc_moves_start",  required_argument, 0,   'm'},
    {"neighbor_dis",  required_argument, 0,   'd'},
    {"initial_tads", no_argument,       0,   'T'},
    {"data_by_input", no_argument,       0,   'D'},
    {"start",  required_argument, 0,   'S'},
    {"end",  required_argument, 0,   'E'},
    {"outfile",  required_argument, 0,   'o'},
    {"out_contact_only",  required_argument, 0,   'C'},
    {"inputfiledir",  required_argument, 0,   'I'},    
    {nullptr, 0, nullptr, 0}
};

static void HelpInfo(char *argv[]) {
    printf("Usage: %s --cfg <cfgfile> or -f <cfgfile>\n", argv[0]);
}

int main(int argc, char *argv[])
{ 
    const char *optstring = "hFDCo:S:E:r:l:t:s:b:m:Td:I:R";   
    int status = 0;
    int arg, index;
    string cfg_file;

    while ((arg = getopt_long(argc, argv, optstring, long_options, &index)) >= 0) {
        if (arg == -1) { 
            HelpInfo(argv); //传递入可执行文件名
            return -1;
        } else {
            // has_options = 1;
            if (arg=='h')
                {HelpInfo(argv);}
            else if (arg=='H'){
                chr=optarg;
            }
            else if (arg=='R'){
                rough_refer=true;
            }            
            else if (arg=='D'){
                data_by_input=true;
            }   
            else if (arg=='S'){
                start_bin=stoi(optarg);
            }   
            else if (arg=='E'){
                end_bin=stoi(optarg);
            }   
            else if (arg=='o'){
                outfile=(optarg);
            }   
            else if (arg=='C'){
                out_contact_only=true;
            }  
            else if (arg=='r'){
                radius=stoi(optarg); 
                assert(radius>0);
            }                                                            
            else if (arg=='t'){
                number_of_threads=stoi(optarg); 
                assert(number_of_threads>0);    
            }              
            else if (arg=='l'){
                learning_rate=stof(optarg); 
                assert(learning_rate>0);
            }    
            else if (arg=='s'){
                update_steps=stoi(optarg);
                assert(update_steps>0);
            }
            else if (arg=='b'){
                burn_in_time=stoi(optarg);
                assert(burn_in_time>0);
            }        
            else if (arg=='m'){
                mc_moves_start=stoi(optarg);
                assert(mc_moves_start>0);
            }     
            else if (arg=='T'){
                initial_tads=true;
            } 
            else if (arg=='d'){
                neighbor_dis=stoi(optarg);
                assert(neighbor_dis>=0);
            }                                         
            else if (arg=='I'){
                inputfiledir=(optarg);
            }  
        }
    }

    reconstruct();
}
// ./test  -m 100000000 -r 40  -d 0  -S 283 -E 12159 -s 5
// ./test --chr 1 -d 5 -r 40 -m 1000000 -o ./ccc.txt -s 4 -I /data/home/txguo/data_use/maxEnt/parameter_access/radius/rough_chr/1/iter_0/m_1000_r_0/polymer