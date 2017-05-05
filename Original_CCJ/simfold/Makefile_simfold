# makefile for the simfold software only

CC = g++
CFLAGS = -g -I./include -I./src/common -I./src/simfold -Wno-deprecated
LFLAGS = -L./ -lsimfold -lm 

SIMFOLD = simfold
SIMFOLD_PF = simfold_pf
TEST_COUNTS = test_get_counts

S_SRCS = simfold.cpp
S_OBJS = $(S_SRCS:.cpp=.o)

S2_SRCS = simfold_pf.cpp
S2_OBJS = $(S2_SRCS:.cpp=.o)

T_SRCS = test_get_counts.cpp
T_OBJS = $(T_SRCS:.cpp=.o)

LIB = libsimfold.a

SRC_C = src/common
SRC_S = src/simfold

L_C_SRCS = $(SRC_C)/params.cpp $(SRC_C)/common.cpp $(SRC_C)/init.cpp $(SRC_C)/timer.cpp

L_S_SRCS = $(SRC_S)/s_specific_functions.cpp $(SRC_S)/s_min_folding.cpp $(SRC_S)/s_hairpin_loop.cpp $(SRC_S)/s_stacked_pair.cpp $(SRC_S)/s_energy_matrix.cpp $(SRC_S)/s_internal_loop.cpp $(SRC_S)/s_multi_loop.cpp $(SRC_S)/s_multi_loop_sub.cpp $(SRC_S)/s_sub_folding.cpp $(SRC_S)/s_partition_function.cpp

L_C_OBJS = $(L_C_SRCS:.cpp=.o)
L_S_OBJS = $(L_S_SRCS:.cpp=.o)

AR = ar
ARFLAGS = rv
RANLIB = ranlib
RM = rm -f

###################################

default: $(LIB) $(SIMFOLD) $(SIMFOLD_PF) $(TEST_COUNTS)

$(LIB): $(L_C_OBJS) $(L_S_OBJS) 
	$(RM) $@
	$(AR) $(ARFLAGS) $@ $(L_C_OBJS) $(L_S_OBJS) 
	$(RANLIB) $@        
   
$(SIMFOLD): $(S_OBJS) $(LIB)
	$(CC) -o $@ $@.o $(LFLAGS)
  
$(SIMFOLD_PF): $(S2_OBJS) $(LIB)
	$(CC) -o $@ $@.o $(LFLAGS)

$(TEST_COUNTS): $(T_OBJS) $(LIB)
	$(CC) -o $@ $@.o $(LFLAGS)
                           
.SUFFIXES: .cpp
.cpp.o:
	$(CC) -c $*.cpp $(CFLAGS) -o $*.o

clean:
	rm -f $(L_C_OBJS) $(L_S_OBJS) 
	rm -f $(S_OBJS) $(S2_OBJS) $(SIMFOLD) $(T_OBJS)*~ core
	rm -f $(SRC_C)/*~ $(SRC_S)/*~ include/*~ *.a *.bak

depend:
	makedepend -- $(CFLAGS) -- $(L_C_SRCS) $(L_S_SRCS) $(S_SRCS) $(S2_SRCS) $(T_SRCS)
	
	
# end of makefile


