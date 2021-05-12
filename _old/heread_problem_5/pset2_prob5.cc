#include "pset2_odes.hh"
#include "pset2_solvers.hh"
#include <iostream>

void version1(double tmax){
    //situation 1:
    agents ag(1250,0.5,0.5);
    //sol_out fsal_solver(double tmax, double lambda, ode *problem, bool writeq, std::string filename = "");
    count_fin pres = fsal_solver(tmax,1e-6,&ag,1,"pset2_prob5_v1.csv");
    std::cout << "function calls: " << pres.first << "\n";
}

void version2(double tmax){
    //situation 1:
    agents ag(1250,0.3,-0.2);
    //sol_out fsal_solver(double tmax, double lambda, ode *problem, bool writeq, std::string filename = "");
    count_fin pres = fsal_solver(tmax,1e-6,&ag,1,"pset2_prob5_v2.csv");
    std::cout << "function calls: " << pres.first << "\n";
}

void version3(double tmax){
    //situation 1:
    agents ag(1250,1,-0.2);
    //sol_out fsal_solver(double tmax, double lambda, ode *problem, bool writeq, std::string filename = "");
    count_fin pres = fsal_solver(tmax,1e-6,&ag,1,"pset2_prob5_v3.csv");
    std::cout << "function calls: " << pres.first << "\n";
}

void testing(double tmax){
    agents ag(1250,0.1,1);
    count_fin pres = fsal_solver(tmax,1e-6,&ag,1,"pset2_prob5_testing.csv");
    std::cout << "function calls: " << pres.first << "\n";
}

void part_a(){
    double tmax = 200;
    //version1(tmax);
    version2(tmax);
    version3(tmax);
    //testing(tmax);
}

void part_b(){
    double tmax = 200;
    agents3d ag3(1000,0.5,0);
    count_fin pres = fsal_solver(tmax,1e-6,&ag3,1,"pset2_prob5_testing.csv");
    std::cout << "function calls: " << pres.first << "\n";
}

int main(){
    part_b();
}
