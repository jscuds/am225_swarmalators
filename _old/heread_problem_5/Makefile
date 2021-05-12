#makefile
#first, specify compiler

cxx = g++
#cflags = -ggdb -fsanitize=address
cflags = 

clean:
	rm pset2_odes.o pset2_solvers.o p1 p7 p5

#build object file for all odes
pset2_odes.o: pset2_odes.cc pset2_odes.hh
	$(cxx) $(cflags) -c pset2_odes.cc

pset2_solvers.o: pset2_solvers.cc pset2_solvers.hh
	$(cxx) $(cflags) -c pset2_solvers.cc

#build prob 1 file
p1: pset2_prob1.cc pset2_odes.o pset2_solvers.o
	$(cxx) $(cflags) -o p1 pset2_prob1.cc pset2_odes.o pset2_solvers.cc

p7: pset2_prob7.cc pset2_odes.o pset2_solvers.o
	$(cxx) $(cflags) -o p7 pset2_prob7.cc pset2_odes.o pset2_solvers.cc

p5: pset2_prob5.cc pset2_odes.o pset2_solvers.o
	$(cxx) $(cflags) -O3 -o p5 pset2_prob5.cc pset2_odes.o pset2_solvers.cc