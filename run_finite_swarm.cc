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
    
    double J = 1, K = 0, N = 1000;
    double r;
    
    for (int i=0; i<5; i++) {
    
        r = 0.5*i;
    
        std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,r));
        o->filename = "figs_data/test_swarm_r"+std::to_string(r)+".txt";
        o->solve_fixed(250.,1e-6,true,500);
        
        std::unique_ptr<swarm_sol_3d> p(new swarm_sol_3d(J,K,N,r));
        p->filename = "figs_data/test_swarm3d_r"+std::to_string(r)+".txt";
        p->solve_fixed(250.,1e-6,true,500);
        
        printf("Finished with r=%g\n",r);
        
    }
    
}

