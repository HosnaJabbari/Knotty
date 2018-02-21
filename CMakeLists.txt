# Ian Wark July 6 2017

cmake_minimum_required(VERSION 3.1)
project(CCJ)

# set SIMFOLD_HOME to simfold folder so CCJ binary can find it
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
set(CMAKE_CXX_FLAGS "-Wall -g -O2")

# SW - deactivate asserts; use for 'production' compile
# add_definitions(-DNDEBUG -Wno-unused-variable -Wno-unused-but-set-variable)

# turn off some warnings
add_definitions(-Wno-format)

add_subdirectory(simfold)
add_subdirectory(Sparse_CCJ)
add_subdirectory(Original_CCJ)
