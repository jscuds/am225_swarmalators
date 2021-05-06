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
    
    double J = 1, K = -0.2, N = 1000;
    
    std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,0.5));
    //swarm_sol* o; o=new swarm_sol(J,K,N,0.5);
    o->filename = "test_finite.txt";
    o->solve_fixed(10,1e-6,true,5);

    //delete o;
    
}

