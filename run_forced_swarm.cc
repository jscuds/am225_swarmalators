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
    
    double J = 1, K = 0, N = 1000, r = -1;
    double F_freq=3*M_PI/2, F_locx=0, F_locy=0, F_locz=0;
    double F;
    
    for (int f=0; f<6; f++) {
        F = double(f);
        
        std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,r,F,F_freq,F_locx,F_locy,F_locz));
        o->filename = "figs_data/forced_swarm_F"+std::to_string(f)+".txt";
        o->solve_fixed(250.,1e-6,true,500);
        
        std::unique_ptr<swarm_sol_3d> p(new swarm_sol_3d(J,K,N,r,F,F_freq,F_locx,F_locy,F_locz));
        p->filename = "figs_data/forced_swarm3d_F"+std::to_string(f)+".txt";
        p->solve_fixed(250.,1e-6,true,500);
        
        printf("Finished with F=%g\n",F);
    }
    
}

