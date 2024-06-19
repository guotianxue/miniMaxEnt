g++ -o test  -I D:\MinGW\include\c++\11.2.0\boost -I D:\MinGW\include\c++\11.2.0\eigen3 .\test.cpp


g++ -std=c++11 -o test  -I /data/software/HiC/vscode_package/boost_1_79_0 -I /data/software/HiC/vscode_package/eigen-3.4.0  ./main.cpp  -pthread
g++ -std=c++11 -o test  ./main.cpp  -pthread

./test

./test  -m 100000000 -r 40  -d 0 -F -S 283 -E 12159 -s 5
#用法见code use python中pbs文件
./test -o /data/home/txguo/data_use/maxEnt/whole_chr/small_tad
# ./test --chr 1 -d 5 -r 40 -m 1000000 -o /data/home/txguo/code/reconstruct_tads_neighbor3  -s 15 -I /data/home/txguo/data_use/maxEnt/whole_chr/chr_tads_structure/1/polymer