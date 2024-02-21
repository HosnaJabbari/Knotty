# Ian Wark March 26 2018

cmake_minimum_required(VERSION 3.1)
project(Knotty)

# set SIMFOLD_HOME to simfold folder so Knotty binary can find it
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)
set(CMAKE_CXX_FLAGS "-Wall -g -O2")

# SW - deactivate asserts; use for 'production' compile
add_definitions(-DNDEBUG -Wno-unused-variable -Wno-unused-but-set-variable)

# turn off some warnings
add_definitions(-Wno-format)

# build simfold
add_subdirectory(simfold)

# build Knotty
SET(SOURCE
  h_externs.h h_globals.h h_struct.h
  simple_map.h
  cmd_line_options.h
  W_final.h W_final.cpp
  pseudo_loop.h pseudo_loop.cpp
  V_final.h V_final.cpp
  VM_final.h VM_final.cpp
  h_common.h h_common.cpp
  trace_arrow.h trace_arrow.cpp
  candidate_list.h candidate_list.cpp
  knotty.h knotty.cpp
)

# pass SIMFOLD_HOME to binary
add_definitions(-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\")

# set path for executable to be created
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executable
add_executable(knotty ${SOURCE})
# set executable name
set_target_properties(knotty PROPERTIES OUTPUT_NAME knotty)

# enable required C++11 features
target_compile_features(knotty PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(knotty LINK_PUBLIC simfold)
