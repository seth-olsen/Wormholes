CC		 = gcc
CXX		 = g++
DEBUG		 = -g
OPTFLAGS	 = -O2
LIBDIR_LAPACK	 = /usr/lib/lapack-3.8.0
LGFORTRAN_PATH	 = /usr/lib/gfortran
INCLUDE		 = -I$(LIBDIR_LAPACK)/LAPACKE/include
CFLAGS		 = -std=c++11 -Wall ${DEBUG} ${OPTFLAGS} $(INCLUDE)
CXXFLAGS         = -std=c++11 -Wall ${DEBUG} ${OPTFLAGS} $(INCLUDE)
LDFLAGS	 	 = -Wall ${DEBUG}
LDLIBS		 = -lbbhutil -L$(LIBDIR_LAPACK) -llapacke -llapack -lblas -L$(LGFORTRAN_PATH) -lgfortran -lm
LOCDIR		 = home/seth/research/simulations

ellis-sim: ellis-sim.o
	-${CXX} -o ellis ellis-sim.o ${LDLIBS}
	rm -f ellis-sim.o

ellis-sim2: ellis-sim.o
	-${CXX} -o ellis2 ellis-sim.o ${LDLIBS}
	rm -f ellis-sim.o

ellis-sim.o: solvers.h ellis-proc.h jacobian.h ellis-fns.h fda-fns.h fda-io.h sim-structs.h sim-header.h sim-init.h ellis-sim.cpp
	$(CXX) -c $(CXXFLAGS) ellis-sim.cpp

ellis-debug: ellis-debug.o
	-${CXX} -o ellis-debug ellis-debug.o ${LDLIBS}
	rm -f ellis-debug.o

ellis-debug.o: solvers.h ellis-proc.h jacobian.h ellis-fns.h fda-fns.h fda-io.h ellis-debug.cpp
	$(CXX) -c $(CXXFLAGS) ellis-debug.cpp

ellis-conv: ellis-conv.o
	$(CXX) -o ellis-conv ellis-conv.o ${LDLIBS}
	rm -f ellis-conv.o

ellis-conv.o: solvers.h ellis-proc.h jacobian.h ellis-fns.h fda-fns.h fda-io.h sim-structs.h sim-header.h sim-init.h ellis-sim.cpp
	$(CXX) -c $(CXXFLAGS) ellis-conv.cpp

ellis-analysis: ellis-analysis.o
	$(CXX) -o ellis-analysis ellis-analysis.o ${LDLIBS}
	rm -f ellis-analysis.o

ellis-analysis.o: solvers.h ellis-proc.h jacobian.h ellis-fns.h fda-fns.h fda-io.h sim-structs.h sim-header.h sim-init.h ellis-sim.cpp
	$(CXX) -c $(CXXFLAGS) ellis-analysis.cpp

ellis-null: ellis-null.o
	$(CXX) -o ellis-null ellis-null.o ${LDLIBS}
	rm -f ellis-null.o

ellis-null.o: solvers.h ellis-proc.h jacobian.h ellis-fns.h fda-fns.h fda-io.h sim-structs.h sim-header.h sim-init.h ellis-sim.cpp
	$(CXX) -c $(CXXFLAGS) ellis-null.cpp

ellis-mass: ellis-mass.o
	$(CXX) -o ellis-mass ellis-mass.o ${LDLIBS}
	rm -f ellis-mass.o

ellis-mass.o: solvers.h ellis-proc.h jacobian.h ellis-fns.h fda-fns.h fda-io.h sim-structs.h sim-header.h sim-init.h ellis-sim.cpp
	$(CXX) -c $(CXXFLAGS) ellis-mass.cpp

test-code: test-code.o
	-${CXX} -o test-code test-code.o ${LDLIBS}
	rm -f test-code.o

test-code.o: ellis-proc.h ellis-fns.h fda-fns.h fda-io.h test-code.cpp
	$(CXX) -c $(CXXFLAGS) test-code.cpp

test-data: test-data.o
	-${CXX} -o test-data test-data.o ${LDLIBS}
	rm -f test-data.o

test-data.o: ellis-proc.h ellis-fns.h fda-fns.h fda-io.h test-data.cpp
	$(CXX) -c $(CXXFLAGS) test-data.cpp

.PHONY : clean
clean :
	rm -f *.o *~

# "make" uses the following implicit rule to compile
# a .c file into a .o file:
#
#       $(CC) -c $(CFLAGS) <the-.c-file>
#
# or to compile a .cpp file into a .o file:
# 
#       $(CXX) -c $(CXXFLAGS) <the-.c-file>
#
# and the following implicit rule for linking:
#
#       $(CC) $(LDFLAGS) <all-dependent-.o-files> $(LDLIBS)
#
# (regardless of whether .o files from C++ or C source) 
# which means the above commands are inserted automatically by
# compiler if no commands given after "target : dependency"
# (note that commands signalled by tab)
#
# file.o implicitly takes file.cpp (or .c) as a dependency,
# so it is inserted by the compiler if not given explicitly
# 
# "make" uses the following rules to decide what to run:
#
#    If you type make and there's a file called makefile,
# it  runs the commands from first target of that file,
# provided dependent files are more recent than target.
#    If you type make and there's a file called Makefile,
# but none called makefile, it runs the commands from first
# target of that file, provided dependent files are more
# recent than target.
#    If you type make -f <filename>  it runs the commands
# from the first target of <filename>, provided dependent files
# are more recent than target.
#    If you type make <target> it looks for a file called
# makefile (then Makefile) and locates the target. 
#
# In all cases it will also run commands for targets that
# are dependencies of the first or given target.
