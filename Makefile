#makefile for final project
cxx = g++
cflags = -Wall -ggdb -fopenmp #-fsanitize=address

execs = test_swarm
objs = fsal_rk4d.o swarm_ode.o

clean:
	rm -f $(execs) $(objs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

test_swarm: test_swarm.cc $(objs)
	$(cxx) $(cflags) -o $@ $^