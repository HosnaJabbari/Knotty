# Ian Wark July 6 2017

cmake_minimum_required(VERSION 2.8)
project(CCJ)

# set SIMFOLD_HOME to simfold folder so CCJ binary can find it
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
add_subdirectory(simfold)
add_subdirectory(Sparse_CCJ)
add_subdirectory(Original_CCJ)
