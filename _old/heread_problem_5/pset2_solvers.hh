#ifndef pset2_solvers
#define pset2_solvers


#include "pset2_odes.hh"
#include <string>

double euclidean_distance(double *y1, double *y2, int n);
double* rk4_final_check(double tmax, int steps, ode *problem);

typedef std::pair<int, double> sol_out;
typedef std::pair<int, double*> count_fin;
count_fin fsal_solver(double tmax, double lambda, ode *problem, bool writeq, std::string filename = "");
sol_out prob1_processing(double tmax, double lambda, ode *problem, bool writeq, std::string filename = "");
//double* fsal_final_check(double tmax, double atol, double rtol, ode *problem);
count_fin fsal_prob7(double tmax, double atol, double rtol, ode *problem);
//void rk4_solver(double tmax, int steps, ode *problem, bool writeq, std::string filename = "");

#endif