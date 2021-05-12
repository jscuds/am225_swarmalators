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

int main() {
    
    double J = 1, K = 0, N = 250, r = -1;
    
    std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,r));
    o->filename = "figs_data/test_swarm.txt";
    o->solve_fixed(250.,1e-6,true,500);
    
    std::unique_ptr<swarm_sol_3d> p(new swarm_sol_3d(J,K,N,r));
    p->filename = "figs_data/test_swarm3d.txt";
    p->solve_fixed(250.,1e-6,true,500);
    
}

