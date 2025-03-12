# Ian Wark July 6 2017

# make sure we set SIMFOLD_HOME correctly to this same directory
# if not, give error
if(${SIMFOLD_HOME} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})

else()
    message(FATAL_ERROR "SIMFOLD_HOME != CMAKE_CURRENT_SOURCE_DIR")
endif()

SET(SOURCE
  src/common/params.cpp 
  src/common/common.cpp
  src/common/init.cpp
  src/simfold/s_specific_functions.cpp 
  src/simfold/s_min_folding.cpp
  src/simfold/s_hairpin_loop.cpp 
  src/simfold/s_stacked_pair.cpp
  src/simfold/s_energy_matrix.cpp 
  src/simfold/s_internal_loop.cpp
  src/simfold/s_multi_loop.cpp 
  src/simfold/s_multi_loop_sub.cpp
  src/simfold/s_sub_folding.cpp 
  src/simfold/s_partition_function.cpp
)
# set location to make .a file
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})

# create library at LIBRARY_OUTPUT_PATH
add_library (simfold ${SOURCE})

# enable required C++11 features
target_compile_features(simfold PRIVATE cxx_nullptr)

# Enable compiler to find files for simfold library
target_include_directories (simfold PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (simfold PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories (simfold PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/common)
target_include_directories (simfold PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src/simfold)

# suppress some warnings
add_definitions(-Wno-format -Wno-unused-result -Wno-sign-compare
  -Wno-maybe-uninitialized -Wno-unused-variable -Wno-unused-but-set-variable)
