# Include the configuration file which defines cxx and cflags
include ./config.mk

# Lists of files to be built
objs= problem_1_rk4d.o
src=$(patsubst %.o,%.cpp,$(objs) $(mg_objs))
execs=p5

all: $(execs)


# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

p5: problem_5.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend