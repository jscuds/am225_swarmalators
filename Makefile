# Makefile for final project
cxx = g++ -fopenmp
cflags = -Wall -ggdb -O3 #-fsanitize=address -ansi -pedantic
clibs = -lgsl

# Lists of files to be built
objs = fsal_rk4d.o swarm_ode.o
src=$(patsubst %.o,%.cc,$(objs) $(mg_objs))
execs = basic_swarm forced_swarm finite_swarm

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(execs) $(objs)

%.o: %.cc
	$(cxx) $(cflags) -c $< $(clibs)

basic_swarm: run_basic_swarm.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(clibs)

forced_swarm: run_forced_swarm.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(clibs)

finite_swarm: run_finite_swarm.cc $(objs)
	$(cxx) $(cflags) -o $@ $^ $(clibs)

.PHONY: clean depend
