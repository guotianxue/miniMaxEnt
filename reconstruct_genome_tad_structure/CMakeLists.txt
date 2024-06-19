cmake_minimum_required(VERSION 3.6)
project(Processed_new_r0)

# set(CMAKE_CXX_STANDARD 11)
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(D:/MinGW/x86_64-w64-mingw32/include/eigen3)
include_directories(D:/MinGW/x86_64-w64-mingw32/include/boost)
# set(SOURCE_FILES main.cpp Initialize.h Moves.h Energy_changes.h global.h Functions.h)
set(SOURCE_FILES test.cpp initial_polymer.h Moves_liner.h  global.h Inverse_Monte_Carlo.h)
# add_definitions("-Wall -g")
add_executable(Raw_data_Inverse ${SOURCE_FILES})