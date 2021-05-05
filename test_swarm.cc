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
    
    double J = 0.5, K = 0.5, N = 1250;
    
    std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N)); //,0.5));
    //swarm_sol* o; o=new swarm_sol(J,K,N,0.5);
    o->filename = "out_J_0.5_K_0.5_var_TEST";
    o->solve_fixed(250.,0.000001,false,500);

    //delete o;
    
}

