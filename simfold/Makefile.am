## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS= subdir-objects

ACLOCAL_AMFLAGS = -I m4

AM_CXXFLAGS += -I$(srcdir)/include -I$(srcdir)/src/simfold -I$(srcdir)/src/common

lib_LTLIBRARIES=libsimfold.la

libsimfold_la_SOURCES = src/common/params.cpp src/common/common.cpp		\
  src/common/init.cpp src/common/timer.cpp src/simfold/s_specific_functions.cpp	\
  src/simfold/s_min_folding.cpp src/simfold/s_hairpin_loop.cpp			\
  src/simfold/s_stacked_pair.cpp src/simfold/s_energy_matrix.cpp		\
  src/simfold/s_internal_loop.cpp src/simfold/s_multi_loop.cpp			\
  src/simfold/s_multi_loop_sub.cpp src/simfold/s_sub_folding.cpp		\
  src/simfold/s_partition_function.cpp src/simfold/s_internal_loop.h		\
  src/simfold/s_min_folding.h src/simfold/s_multi_loop_sub.h			\
  src/simfold/s_partition_function_RESTRICTED.h				\
  src/simfold/s_partition_function_UNRESTRICTED.h				\
  src/simfold/s_specific_functions.h					\
  src/simfold/s_partition_function_complex.h src/simfold/s_stacked_pair.h	\
  src/simfold/s_hairpin_loop.h src/simfold/s_partition_function.h		\
  src/simfold/s_multi_loop.h src/simfold/s_sub_folding.h			\
  src/simfold/s_energy_matrix.h src/common/timer.h src/common/common.h		\
  src/common/params.h


# link binaries to libsimfold
LDADD = libsimfold.la

bin_PROGRAMS = simfold simfold_pf test_get_counts

simfold_SOURCES=simfold.cpp
simfold_pf_SOURCES=simfold_pf.cpp
test_get_counts_SOURCES=test_get_counts.cpp
