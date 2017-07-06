# Ian Wark July 6 2017

SET(SOURCE
  h_externs.h h_globals.h h_struct.h	
  base.h simple_map.h 
  cmd_line_options.h cmd_line_options.cpp  
  W_final.h W_final.cpp 
  pseudo_loop.h pseudo_loop.cpp 
  V_final.h V_final.cpp 
  VM_final.h VM_final.cpp
  h_common.h h_common.cpp 
  trace_arrow.h trace_arrow.cpp  	
  candidate_list.h candidate_list.cpp 
  CCJ.h CCJ.cpp 
)

# enable c++11
add_compile_options(-std=c++11)
# pass SIMFOLD_HOME to binary
add_definitions("-DSIMFOLD_HOME=\"${SIMFOLD_HOME}\"")

# set path for executable to be created
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
# create executable
add_executable(Sparse_CCJ ${SOURCE})
set_target_properties(Sparse_CCJ PROPERTIES OUTPUT_NAME CCJ)

# link to simfold
target_link_libraries(Sparse_CCJ LINK_PUBLIC simfold)
