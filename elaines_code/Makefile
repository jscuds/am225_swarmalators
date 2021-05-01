# Include the configuration file

# Lists of files to be built
objs= problem_1_rk4d.o
src=$(patsubst %.o,%.cpp,$(objs) $(mg_objs))
execs=p5

all: $(execs)


# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	g++-10 -fopenmp -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cpp
	g++-10 -g -Wall -Wextra -fopenmp -c $<

p5: problem_5.cpp $(objs)
	g++-10 -g -Wall -Wextra -fopenmp -O3 -o p5 problem_5.cpp problem_1_rk4d.o 

.PHONY: clean depend
