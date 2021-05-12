# AM225 Final Project (Spring 2021)
# Exploration of Efficient Simulation for Swarmalators
**Contributors: M. Elaine Cunha, Helen Read, Jack Scudder**

See our report in `final_paper/` for more information. 

## Directory
- `_old/`:  folder containing out-of-date code and figures
- `figs/`:  folder with figures used in report
- `figs_data/`:  folder with data for producing figures in report (incomplete as some files are quite large)
- `final_paper/`:  folder containing our final project report
- `Figures.ipynb`:  Jupyter notebook for producing figures in report
- `Makefile`:  compiles code; check that line 2 contains the correct compiler for your system before running; type `make` in the command line to compile all code
- `Makefile.dep`:  details dependencies for compiling code
- `fsal_rk4d.cc`:  implements five-step Runge-Kutta first-same-as-last adaptive integration with dense output
- `fsal_rk4d.hh`:  header file for `fsal_rk4d.cc`
- `point.hh`:  header file for a C++ struct used for storing agent data in the grid-binning implementation of our finite-cutoff scenario
- `run_basic_swarm.cc`:  outputs 2D and 3D swarmalator results for our base case (after typing `make`, use `./basic_swarm` to run); see file for swarmalator parameters and assumptions
- `run_finite_swarm.cc`:  outputs 2D and 3D swarmalator results with a finite-cutoff interaction radius (after typing `make`, use `./finite_swarm` to run); see file for swarmalator parameters and assumptions
- `run_forced_swarm.cc`:  outputs 2D and 3D swarmalator results with external forcing (after typing `make`, use `./forced_swarm` to run); see file for swarmalator parameters and assumptions
- `swarm_ode.cc`:  implements swarmalator differential equation with variations
- `swarm_ode.hh`:  header file for `swarm_ode.cc`
