ACLOCAL_AMFLAGS  = -I m4

bin_PROGRAMS=CCJ

test = -Wl,-rpath,$(LIBPATH)
CCJ_LDFLAGS = $(test)
CCJ_LDADD =  $(SIMFOLD_LIBS) 

CCJ_SOURCES = \
	h_externs.h h_globals.h h_struct.h	\
	base.h simple_map.h \
	cmd_line_options.h cmd_line_options.cpp  \
	W_final.h W_final.cpp \
	pseudo_loop.h pseudo_loop.cpp \
	V_final.h V_final.cpp \
	VM_final.h VM_final.cpp \
	h_common.h h_common.cpp \
	trace_arrow.h trace_arrow.cpp  	\
	candidate_list.h candidate_list.cpp \
	CCJ.h CCJ.cpp 
