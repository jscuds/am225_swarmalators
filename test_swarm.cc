//
//  problem_5.cpp
//  hw2_code
//
//  Created by Elaine Cunha on 3/3/21.
//
#include <memory>

#include "swarm_ode.hh"
#include "fsal_rk4d.hh"

/** The number of components in the ODE system. */
//const int ns=3;

int main() {
    
    double J = 0.5, K = 0.5, N = 250, r = -1.;
//    int dim = 2;
    
    std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,r));
    //swarm_sol* o; o=new swarm_sol(J,K,N,0.5);
    o->filename = "test_finite.txt";
    o->solve_fixed(250.,1e-6,true,500);
    
}

