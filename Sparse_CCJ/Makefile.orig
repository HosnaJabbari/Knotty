# makefile for the HFold software only

CC = g++
CFLAGS = -g  -I./simfold/include -I./simfold/src/common -I./simfold/src/simfold -I ./ -Wno-deprecated -std=c++11 -DNDEBUG
LFLAGS = -L./ -lsimfold -lm

SIMFOLD = simfold
SIMFOLD_PF = simfold_pf
CCJ = CCJ


S_SRCS = simfold.cpp
S_OBJS = $(S_SRCS:.cpp=.o)

S2_SRCS = simfold_pf.cpp
S_OBJS = $(S_SRCS:.cpp=.o)

C_SRCS = CCJ.cpp
C_OBJS = $(C_SRCS:.cpp=.o)

LIB = libsimfold.a

SRC_C = simfold/src/common
SRC_S = simfold/src/simfold

L_C_SRCS = $(SRC_C)/params.cpp $(SRC_C)/common.cpp $(SRC_C)/init.cpp $(SRC_C)/timer.cpp

L_S_SRCS = $(SRC_S)/s_specific_functions.cpp $(SRC_S)/s_min_folding.cpp $(SRC_S)/s_hairpin_loop.cpp $(SRC_S)/s_stacked_pair.cpp $(SRC_S)/s_energy_matrix.cpp $(SRC_S)/s_internal_loop.cpp $(SRC_S)/s_multi_loop.cpp $(SRC_S)/s_multi_loop_sub.cpp $(SRC_S)/s_sub_folding.cpp $(SRC_S)/s_partition_function.cpp 

L_H_SRCS = cmd_line_options.cpp W_final.cpp pseudo_loop.cpp V_final.cpp VM_final.cpp h_common.cpp trace_arrow.cpp candidate_list.cpp

L_C_OBJS = $(L_C_SRCS:.cpp=.o)
L_S_OBJS = $(L_S_SRCS:.cpp=.o)
L_H_OBJS = $(L_H_SRCS:.cpp=.o)

AR = ar
ARFLAGS = rv
RANLIB = ranlib
RM = rm -f

###################################

default: $(LIB) $(CCJ) 

$(LIB): $(L_C_OBJS) $(L_S_OBJS) $(L_H_OBJS) 
	$(RM) $@
	$(AR) $(ARFLAGS) $@ $(L_C_OBJS) $(L_S_OBJS) $(L_H_OBJS)
	$(RANLIB) $@        
   
$(SIMFOLD): $(S_OBJS) $(LIB)
	$(CC) -o $@ $@.o $(LFLAGS)
                             
$(SIMFOLD_PF): $(S2_OBJS) $(LIB)
	$(CC) -o $@ $@.o $(LFLAGS)

$(CCJ): $(C_OBJS) $(LIB)
	$(CC) -o $@ $@.o $(LFLAGS)

	
.SUFFIXES: .cpp
.cpp.o:
	$(CC) -c $*.cpp $(CFLAGS) -o $*.o

clean:
	rm -f $(L_C_OBJS) $(L_S_OBJS) $(L_H_OBJS)
	#rm -f $(S_OBJS) $(SIMFOLD) *~ core
	rm -f $(C_OBJS) $(CCJ) *~ core
	#rm -f $(CR_OBJS) *~ core
	#rm -f $(PA_OBJS) *~ core
	rm -f $(SRC_C)/*~ $(SRC_S)/*~ include/*~ *.a *.bak

depend:
	makedepend -- $(CFLAGS) -- $(L_C_SRCS) $(L_S_SRCS) $(S_SRCS) $(S2_SRCS) 
	
	
# end of makefile


# DO NOT DELETE

simfold/src/common/params.o: /usr/include/stdio.h /usr/include/features.h
simfold/src/common/params.o: /usr/include/stdc-predef.h
simfold/src/common/params.o: /usr/include/sys/cdefs.h
simfold/src/common/params.o: /usr/include/bits/wordsize.h
simfold/src/common/params.o: /usr/include/gnu/stubs.h
simfold/src/common/params.o: /usr/include/bits/types.h
simfold/src/common/params.o: /usr/include/bits/typesizes.h
simfold/src/common/params.o: /usr/include/libio.h /usr/include/_G_config.h
simfold/src/common/params.o: /usr/include/wchar.h
simfold/src/common/params.o: /usr/include/bits/stdio_lim.h
simfold/src/common/params.o: /usr/include/bits/sys_errlist.h
simfold/src/common/params.o: /usr/include/math.h /usr/include/bits/huge_val.h
simfold/src/common/params.o: /usr/include/bits/huge_valf.h
simfold/src/common/params.o: /usr/include/bits/huge_vall.h
simfold/src/common/params.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
simfold/src/common/params.o: /usr/include/bits/mathdef.h
simfold/src/common/params.o: /usr/include/bits/mathcalls.h
simfold/src/common/params.o: /usr/include/stdlib.h
simfold/src/common/params.o: /usr/include/bits/waitflags.h
simfold/src/common/params.o: /usr/include/bits/waitstatus.h
simfold/src/common/params.o: /usr/include/endian.h /usr/include/bits/endian.h
simfold/src/common/params.o: /usr/include/bits/byteswap.h
simfold/src/common/params.o: /usr/include/bits/byteswap-16.h
simfold/src/common/params.o: /usr/include/sys/types.h /usr/include/time.h
simfold/src/common/params.o: /usr/include/sys/select.h
simfold/src/common/params.o: /usr/include/bits/select.h
simfold/src/common/params.o: /usr/include/bits/sigset.h
simfold/src/common/params.o: /usr/include/bits/time.h
simfold/src/common/params.o: /usr/include/sys/sysmacros.h
simfold/src/common/params.o: /usr/include/bits/pthreadtypes.h
simfold/src/common/params.o: /usr/include/alloca.h
simfold/src/common/params.o: /usr/include/bits/stdlib-float.h
simfold/src/common/params.o: /usr/include/ctype.h /usr/include/xlocale.h
simfold/src/common/params.o: simfold/src/common/params.h
simfold/src/common/params.o: ./simfold/include/constants.h
simfold/src/common/params.o: simfold/src/common/common.h
simfold/src/common/params.o: ./simfold/include/structs.h
simfold/src/common/params.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/common/params.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/common/params.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/common/params.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/common/params.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/common/params.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/common/params.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/common/params.o: ./simfold/include/externs.h
simfold/src/common/params.o: /usr/include/string.h
simfold/src/common/params.o: ./simfold/src/simfold/s_specific_functions.h
simfold/src/common/params.o: ./simfold/src/simfold/s_min_folding.h
simfold/src/common/params.o: ./simfold/src/simfold/s_sub_folding.h
simfold/src/common/params.o: ./simfold/include/simfold.h
simfold/src/common/params.o: ./simfold/include/init.h
simfold/src/common/common.o: /usr/include/stdio.h /usr/include/features.h
simfold/src/common/common.o: /usr/include/stdc-predef.h
simfold/src/common/common.o: /usr/include/sys/cdefs.h
simfold/src/common/common.o: /usr/include/bits/wordsize.h
simfold/src/common/common.o: /usr/include/gnu/stubs.h
simfold/src/common/common.o: /usr/include/bits/types.h
simfold/src/common/common.o: /usr/include/bits/typesizes.h
simfold/src/common/common.o: /usr/include/libio.h /usr/include/_G_config.h
simfold/src/common/common.o: /usr/include/wchar.h
simfold/src/common/common.o: /usr/include/bits/stdio_lim.h
simfold/src/common/common.o: /usr/include/bits/sys_errlist.h
simfold/src/common/common.o: /usr/include/stdlib.h
simfold/src/common/common.o: /usr/include/bits/waitflags.h
simfold/src/common/common.o: /usr/include/bits/waitstatus.h
simfold/src/common/common.o: /usr/include/endian.h /usr/include/bits/endian.h
simfold/src/common/common.o: /usr/include/bits/byteswap.h
simfold/src/common/common.o: /usr/include/bits/byteswap-16.h
simfold/src/common/common.o: /usr/include/sys/types.h /usr/include/time.h
simfold/src/common/common.o: /usr/include/sys/select.h
simfold/src/common/common.o: /usr/include/bits/select.h
simfold/src/common/common.o: /usr/include/bits/sigset.h
simfold/src/common/common.o: /usr/include/bits/time.h
simfold/src/common/common.o: /usr/include/sys/sysmacros.h
simfold/src/common/common.o: /usr/include/bits/pthreadtypes.h
simfold/src/common/common.o: /usr/include/alloca.h
simfold/src/common/common.o: /usr/include/bits/stdlib-float.h
simfold/src/common/common.o: /usr/include/math.h /usr/include/bits/huge_val.h
simfold/src/common/common.o: /usr/include/bits/huge_valf.h
simfold/src/common/common.o: /usr/include/bits/huge_vall.h
simfold/src/common/common.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
simfold/src/common/common.o: /usr/include/bits/mathdef.h
simfold/src/common/common.o: /usr/include/bits/mathcalls.h
simfold/src/common/common.o: /usr/include/string.h /usr/include/xlocale.h
simfold/src/common/common.o: ./simfold/include/constants.h
simfold/src/common/common.o: ./simfold/include/structs.h
simfold/src/common/common.o: ./simfold/include/externs.h
simfold/src/common/common.o: simfold/src/common/common.h
simfold/src/common/common.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/common/common.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/common/common.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/common/common.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/common/common.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/common/common.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/common/common.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/common/init.o: /usr/include/stdio.h /usr/include/features.h
simfold/src/common/init.o: /usr/include/stdc-predef.h
simfold/src/common/init.o: /usr/include/sys/cdefs.h
simfold/src/common/init.o: /usr/include/bits/wordsize.h
simfold/src/common/init.o: /usr/include/gnu/stubs.h /usr/include/bits/types.h
simfold/src/common/init.o: /usr/include/bits/typesizes.h /usr/include/libio.h
simfold/src/common/init.o: /usr/include/_G_config.h /usr/include/wchar.h
simfold/src/common/init.o: /usr/include/bits/stdio_lim.h
simfold/src/common/init.o: /usr/include/bits/sys_errlist.h
simfold/src/common/init.o: /usr/include/stdlib.h
simfold/src/common/init.o: /usr/include/bits/waitflags.h
simfold/src/common/init.o: /usr/include/bits/waitstatus.h
simfold/src/common/init.o: /usr/include/endian.h /usr/include/bits/endian.h
simfold/src/common/init.o: /usr/include/bits/byteswap.h
simfold/src/common/init.o: /usr/include/bits/byteswap-16.h
simfold/src/common/init.o: /usr/include/sys/types.h /usr/include/time.h
simfold/src/common/init.o: /usr/include/sys/select.h
simfold/src/common/init.o: /usr/include/bits/select.h
simfold/src/common/init.o: /usr/include/bits/sigset.h
simfold/src/common/init.o: /usr/include/bits/time.h
simfold/src/common/init.o: /usr/include/sys/sysmacros.h
simfold/src/common/init.o: /usr/include/bits/pthreadtypes.h
simfold/src/common/init.o: /usr/include/alloca.h
simfold/src/common/init.o: /usr/include/bits/stdlib-float.h
simfold/src/common/init.o: /usr/include/string.h /usr/include/xlocale.h
simfold/src/common/init.o: /usr/include/math.h /usr/include/bits/huge_val.h
simfold/src/common/init.o: /usr/include/bits/huge_valf.h
simfold/src/common/init.o: /usr/include/bits/huge_vall.h
simfold/src/common/init.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
simfold/src/common/init.o: /usr/include/bits/mathdef.h
simfold/src/common/init.o: /usr/include/bits/mathcalls.h /usr/include/ctype.h
simfold/src/common/init.o: ./simfold/include/constants.h
simfold/src/common/init.o: ./simfold/include/structs.h
simfold/src/common/init.o: ./simfold/include/globals.h
simfold/src/common/init.o: simfold/src/common/common.h
simfold/src/common/init.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/common/init.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/common/init.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/common/init.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/common/init.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/common/init.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/common/init.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/common/init.o: simfold/src/common/params.h
simfold/src/common/timer.o: /usr/include/math.h /usr/include/features.h
simfold/src/common/timer.o: /usr/include/stdc-predef.h
simfold/src/common/timer.o: /usr/include/sys/cdefs.h
simfold/src/common/timer.o: /usr/include/bits/wordsize.h
simfold/src/common/timer.o: /usr/include/gnu/stubs.h
simfold/src/common/timer.o: /usr/include/bits/huge_val.h
simfold/src/common/timer.o: /usr/include/bits/huge_valf.h
simfold/src/common/timer.o: /usr/include/bits/huge_vall.h
simfold/src/common/timer.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
simfold/src/common/timer.o: /usr/include/bits/mathdef.h
simfold/src/common/timer.o: /usr/include/bits/mathcalls.h
simfold/src/common/timer.o: /usr/include/ctype.h /usr/include/bits/types.h
simfold/src/common/timer.o: /usr/include/bits/typesizes.h
simfold/src/common/timer.o: /usr/include/endian.h /usr/include/bits/endian.h
simfold/src/common/timer.o: /usr/include/bits/byteswap.h
simfold/src/common/timer.o: /usr/include/bits/byteswap-16.h
simfold/src/common/timer.o: /usr/include/xlocale.h /usr/include/time.h
simfold/src/common/timer.o: /usr/include/string.h /usr/include/sys/times.h
simfold/src/common/timer.o: /usr/include/limits.h
simfold/src/common/timer.o: /usr/include/bits/posix1_lim.h
simfold/src/common/timer.o: /usr/include/bits/local_lim.h
simfold/src/common/timer.o: /usr/include/linux/limits.h
simfold/src/common/timer.o: /usr/include/bits/posix2_lim.h
simfold/src/common/timer.o: simfold/src/common/timer.h
simfold/src/simfold/s_specific_functions.o: /usr/include/stdio.h
simfold/src/simfold/s_specific_functions.o: /usr/include/features.h
simfold/src/simfold/s_specific_functions.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_specific_functions.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_specific_functions.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/types.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_specific_functions.o: /usr/include/libio.h
simfold/src/simfold/s_specific_functions.o: /usr/include/_G_config.h
simfold/src/simfold/s_specific_functions.o: /usr/include/wchar.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_specific_functions.o: /usr/include/stdlib.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/waitflags.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/waitstatus.h
simfold/src/simfold/s_specific_functions.o: /usr/include/endian.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/endian.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/byteswap.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/byteswap-16.h
simfold/src/simfold/s_specific_functions.o: /usr/include/sys/types.h
simfold/src/simfold/s_specific_functions.o: /usr/include/time.h
simfold/src/simfold/s_specific_functions.o: /usr/include/sys/select.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/select.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/sigset.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/time.h
simfold/src/simfold/s_specific_functions.o: /usr/include/sys/sysmacros.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/pthreadtypes.h
simfold/src/simfold/s_specific_functions.o: /usr/include/alloca.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/stdlib-float.h
simfold/src/simfold/s_specific_functions.o: /usr/include/math.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/huge_val.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/huge_valf.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/huge_vall.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/inf.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/nan.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/mathdef.h
simfold/src/simfold/s_specific_functions.o: /usr/include/bits/mathcalls.h
simfold/src/simfold/s_specific_functions.o: /usr/include/string.h
simfold/src/simfold/s_specific_functions.o: /usr/include/xlocale.h
simfold/src/simfold/s_specific_functions.o: ./simfold/include/constants.h
simfold/src/simfold/s_specific_functions.o: ./simfold/include/structs.h
simfold/src/simfold/s_specific_functions.o: ./simfold/include/externs.h
simfold/src/simfold/s_specific_functions.o: simfold/src/common/common.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_specific_functions.o: ./simfold/include/simfold.h
simfold/src/simfold/s_specific_functions.o: ./simfold/include/init.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_specific_functions.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_min_folding.h
simfold/src/simfold/s_specific_functions.o: ./simfold/src/simfold/s_sub_folding.h
simfold/src/simfold/s_min_folding.o: /usr/include/stdio.h
simfold/src/simfold/s_min_folding.o: /usr/include/features.h
simfold/src/simfold/s_min_folding.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_min_folding.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_min_folding.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/types.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_min_folding.o: /usr/include/libio.h
simfold/src/simfold/s_min_folding.o: /usr/include/_G_config.h
simfold/src/simfold/s_min_folding.o: /usr/include/wchar.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_min_folding.o: /usr/include/string.h
simfold/src/simfold/s_min_folding.o: /usr/include/xlocale.h
simfold/src/simfold/s_min_folding.o: /usr/include/stdlib.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/waitflags.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/waitstatus.h
simfold/src/simfold/s_min_folding.o: /usr/include/endian.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/endian.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/byteswap.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/byteswap-16.h
simfold/src/simfold/s_min_folding.o: /usr/include/sys/types.h
simfold/src/simfold/s_min_folding.o: /usr/include/time.h
simfold/src/simfold/s_min_folding.o: /usr/include/sys/select.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/select.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/sigset.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/time.h
simfold/src/simfold/s_min_folding.o: /usr/include/sys/sysmacros.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/pthreadtypes.h
simfold/src/simfold/s_min_folding.o: /usr/include/alloca.h
simfold/src/simfold/s_min_folding.o: /usr/include/bits/stdlib-float.h
simfold/src/simfold/s_min_folding.o: ./simfold/include/constants.h
simfold/src/simfold/s_min_folding.o: ./simfold/include/structs.h
simfold/src/simfold/s_min_folding.o: ./simfold/include/externs.h
simfold/src/simfold/s_min_folding.o: simfold/src/common/common.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_min_folding.o: ./simfold/include/simfold.h
simfold/src/simfold/s_min_folding.o: ./simfold/include/init.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_specific_functions.h
simfold/src/simfold/s_min_folding.o: ./simfold/src/simfold/s_min_folding.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/string.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/features.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/xlocale.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/ctype.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/types.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/endian.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/endian.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/byteswap.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/byteswap-16.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/math.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/huge_val.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/huge_valf.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/huge_vall.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/inf.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/nan.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/mathdef.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/mathcalls.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/stdlib.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/waitflags.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/waitstatus.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/sys/types.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/time.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/sys/select.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/select.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/sigset.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/time.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/sys/sysmacros.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/pthreadtypes.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/alloca.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/stdlib-float.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/include/constants.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/include/structs.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/include/externs.h
simfold/src/simfold/s_hairpin_loop.o: simfold/src/common/common.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/stdio.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/libio.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/_G_config.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/wchar.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_hairpin_loop.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/include/simfold.h
simfold/src/simfold/s_hairpin_loop.o: ./simfold/include/init.h
simfold/src/simfold/s_hairpin_loop.o: simfold/src/common/params.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/string.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/features.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/xlocale.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/include/constants.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/include/structs.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/include/externs.h
simfold/src/simfold/s_stacked_pair.o: simfold/src/common/common.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/stdio.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/bits/types.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/libio.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/_G_config.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/wchar.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_stacked_pair.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/include/simfold.h
simfold/src/simfold/s_stacked_pair.o: ./simfold/include/init.h
simfold/src/simfold/s_stacked_pair.o: simfold/src/common/params.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/stdio.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/features.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/types.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/libio.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/_G_config.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/wchar.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/string.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/xlocale.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/stdlib.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/waitflags.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/waitstatus.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/endian.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/endian.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/byteswap.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/byteswap-16.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/sys/types.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/time.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/sys/select.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/select.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/sigset.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/time.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/sys/sysmacros.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/pthreadtypes.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/alloca.h
simfold/src/simfold/s_energy_matrix.o: /usr/include/bits/stdlib-float.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/include/constants.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/include/structs.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/include/externs.h
simfold/src/simfold/s_energy_matrix.o: simfold/src/common/common.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/include/simfold.h
simfold/src/simfold/s_energy_matrix.o: ./simfold/include/init.h
simfold/src/simfold/s_internal_loop.o: /usr/include/stdio.h
simfold/src/simfold/s_internal_loop.o: /usr/include/features.h
simfold/src/simfold/s_internal_loop.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_internal_loop.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_internal_loop.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/types.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_internal_loop.o: /usr/include/libio.h
simfold/src/simfold/s_internal_loop.o: /usr/include/_G_config.h
simfold/src/simfold/s_internal_loop.o: /usr/include/wchar.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_internal_loop.o: /usr/include/string.h
simfold/src/simfold/s_internal_loop.o: /usr/include/xlocale.h
simfold/src/simfold/s_internal_loop.o: /usr/include/stdlib.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/waitflags.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/waitstatus.h
simfold/src/simfold/s_internal_loop.o: /usr/include/endian.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/endian.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/byteswap.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/byteswap-16.h
simfold/src/simfold/s_internal_loop.o: /usr/include/sys/types.h
simfold/src/simfold/s_internal_loop.o: /usr/include/time.h
simfold/src/simfold/s_internal_loop.o: /usr/include/sys/select.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/select.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/sigset.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/time.h
simfold/src/simfold/s_internal_loop.o: /usr/include/sys/sysmacros.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/pthreadtypes.h
simfold/src/simfold/s_internal_loop.o: /usr/include/alloca.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/stdlib-float.h
simfold/src/simfold/s_internal_loop.o: /usr/include/ctype.h
simfold/src/simfold/s_internal_loop.o: /usr/include/math.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/huge_val.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/huge_valf.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/huge_vall.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/inf.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/nan.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/mathdef.h
simfold/src/simfold/s_internal_loop.o: /usr/include/bits/mathcalls.h
simfold/src/simfold/s_internal_loop.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_internal_loop.o: ./simfold/include/constants.h
simfold/src/simfold/s_internal_loop.o: ./simfold/include/structs.h
simfold/src/simfold/s_internal_loop.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_internal_loop.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_internal_loop.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_internal_loop.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_internal_loop.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_internal_loop.o: ./simfold/include/externs.h
simfold/src/simfold/s_internal_loop.o: simfold/src/common/common.h
simfold/src/simfold/s_internal_loop.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_internal_loop.o: ./simfold/include/simfold.h
simfold/src/simfold/s_internal_loop.o: ./simfold/include/init.h
simfold/src/simfold/s_internal_loop.o: simfold/src/common/params.h
simfold/src/simfold/s_multi_loop.o: /usr/include/stdio.h
simfold/src/simfold/s_multi_loop.o: /usr/include/features.h
simfold/src/simfold/s_multi_loop.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_multi_loop.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_multi_loop.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_multi_loop.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_multi_loop.o: /usr/include/bits/types.h
simfold/src/simfold/s_multi_loop.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_multi_loop.o: /usr/include/libio.h
simfold/src/simfold/s_multi_loop.o: /usr/include/_G_config.h
simfold/src/simfold/s_multi_loop.o: /usr/include/wchar.h
simfold/src/simfold/s_multi_loop.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_multi_loop.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_multi_loop.o: ./simfold/include/externs.h
simfold/src/simfold/s_multi_loop.o: ./simfold/include/structs.h
simfold/src/simfold/s_multi_loop.o: ./simfold/include/constants.h
simfold/src/simfold/s_multi_loop.o: simfold/src/common/common.h
simfold/src/simfold/s_multi_loop.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_multi_loop.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_multi_loop.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_multi_loop.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_multi_loop.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_multi_loop.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_multi_loop.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_multi_loop.o: ./simfold/include/simfold.h
simfold/src/simfold/s_multi_loop.o: ./simfold/include/init.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/stdio.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/features.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/bits/types.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/libio.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/_G_config.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/wchar.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_multi_loop_sub.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/include/externs.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/include/structs.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/include/constants.h
simfold/src/simfold/s_multi_loop_sub.o: simfold/src/common/common.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_multi_loop_sub.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_sub_folding.o: ./simfold/src/simfold/s_sub_folding.h
simfold/src/simfold/s_sub_folding.o: /usr/include/stdio.h
simfold/src/simfold/s_sub_folding.o: /usr/include/features.h
simfold/src/simfold/s_sub_folding.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_sub_folding.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_sub_folding.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_sub_folding.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_sub_folding.o: /usr/include/bits/types.h
simfold/src/simfold/s_sub_folding.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_sub_folding.o: /usr/include/libio.h
simfold/src/simfold/s_sub_folding.o: /usr/include/_G_config.h
simfold/src/simfold/s_sub_folding.o: /usr/include/wchar.h
simfold/src/simfold/s_sub_folding.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_sub_folding.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_sub_folding.o: /usr/include/string.h
simfold/src/simfold/s_sub_folding.o: /usr/include/xlocale.h
simfold/src/simfold/s_sub_folding.o: ./simfold/include/structs.h
simfold/src/simfold/s_sub_folding.o: ./simfold/include/constants.h
simfold/src/simfold/s_sub_folding.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_sub_folding.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_sub_folding.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_sub_folding.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_sub_folding.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_sub_folding.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_sub_folding.o: ./simfold/include/externs.h
simfold/src/simfold/s_sub_folding.o: simfold/src/common/common.h
simfold/src/simfold/s_sub_folding.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_sub_folding.o: ./simfold/include/simfold.h
simfold/src/simfold/s_sub_folding.o: ./simfold/include/init.h
simfold/src/simfold/s_partition_function.o: /usr/include/stdio.h
simfold/src/simfold/s_partition_function.o: /usr/include/features.h
simfold/src/simfold/s_partition_function.o: /usr/include/stdc-predef.h
simfold/src/simfold/s_partition_function.o: /usr/include/sys/cdefs.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/wordsize.h
simfold/src/simfold/s_partition_function.o: /usr/include/gnu/stubs.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/types.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/typesizes.h
simfold/src/simfold/s_partition_function.o: /usr/include/libio.h
simfold/src/simfold/s_partition_function.o: /usr/include/_G_config.h
simfold/src/simfold/s_partition_function.o: /usr/include/wchar.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/stdio_lim.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/sys_errlist.h
simfold/src/simfold/s_partition_function.o: /usr/include/string.h
simfold/src/simfold/s_partition_function.o: /usr/include/xlocale.h
simfold/src/simfold/s_partition_function.o: /usr/include/stdlib.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/waitflags.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/waitstatus.h
simfold/src/simfold/s_partition_function.o: /usr/include/endian.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/endian.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/byteswap.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/byteswap-16.h
simfold/src/simfold/s_partition_function.o: /usr/include/sys/types.h
simfold/src/simfold/s_partition_function.o: /usr/include/time.h
simfold/src/simfold/s_partition_function.o: /usr/include/sys/select.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/select.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/sigset.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/time.h
simfold/src/simfold/s_partition_function.o: /usr/include/sys/sysmacros.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/pthreadtypes.h
simfold/src/simfold/s_partition_function.o: /usr/include/alloca.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/stdlib-float.h
simfold/src/simfold/s_partition_function.o: /usr/include/math.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/huge_val.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/huge_valf.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/huge_vall.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/inf.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/nan.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/mathdef.h
simfold/src/simfold/s_partition_function.o: /usr/include/bits/mathcalls.h
simfold/src/simfold/s_partition_function.o: /usr/include/ctype.h
simfold/src/simfold/s_partition_function.o: ./simfold/include/constants.h
simfold/src/simfold/s_partition_function.o: ./simfold/include/structs.h
simfold/src/simfold/s_partition_function.o: ./simfold/include/externs.h
simfold/src/simfold/s_partition_function.o: simfold/src/common/common.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_partition_function.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_hairpin_loop.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_internal_loop.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_energy_matrix.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_stacked_pair.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_multi_loop.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_multi_loop_sub.h
simfold/src/simfold/s_partition_function.o: ./simfold/include/simfold.h
simfold/src/simfold/s_partition_function.o: ./simfold/include/init.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_min_folding.h
simfold/src/simfold/s_partition_function.o: ./simfold/src/simfold/s_sub_folding.h
simfold/src/simfold/s_partition_function.o: simfold/src/common/params.h
