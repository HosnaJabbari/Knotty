# Ian Wark March 26 2018

cmake_minimum_required(VERSION 3.5)
project(Knotty)

# set SIMFOLD_HOME to simfold folder so Knotty binary can find it
set(SIMFOLD_HOME ${CMAKE_SOURCE_DIR}/simfold)


# debugging flags (makes the program significantly slower)
# -fsanitize=address    # out-of-bounds reads/writes
# -fsanitize=undefined  # undefined behavior (e.g. int overflow)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O3")


# SW - deactivate asserts; use for 'production' compile
add_definitions(-DNDEBUG -Wno-unused-variable -Wno-unused-but-set-variable)

# turn off some warnings
add_definitions(-Wno-format)

# build simfold
add_subdirectory(simfold)

# build Knotty
SET(SOURCE
  src/W_final.cpp
  src/pseudo_loop.cpp
  src/V_final.cpp
  src/VM_final.cpp
  src/h_common.cpp
  src/trace_arrow.cpp
  src/candidate_list.cpp
  src/knotty.cpp
  src/param_path.cpp
)

# pass SIMFOLD_HOME to binary
add_definitions(-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\")

# set path for executable to be created
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executable
add_executable(knotty ${SOURCE})
# set executable name
set_target_properties(knotty PROPERTIES OUTPUT_NAME knotty)

# Install Knotty into the bin directory (Needed for conda)
install(TARGETS knotty RUNTIME DESTINATION bin)

# enable required C++11 features
target_compile_features(knotty PRIVATE cxx_nullptr)

# link to simfold
target_link_libraries(knotty LINK_PUBLIC simfold)
