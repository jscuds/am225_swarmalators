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
    
    double J = 0.1, K = 1, N = 100;
    double r[4];
    r[0]=-1; r[1]=2; r[2]=1.;r[3]=0.5;
    
    for (int i=0; i<4; i++) {
        
        std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,r[i]));
        o->filename = "figs_data/finite_swarm_r"+std::to_string(r[i])+".txt";
        o->solve_fixed(250.,1e-6,true,500);
        
        std::unique_ptr<swarm_sol_3d> p(new swarm_sol_3d(J,K,N,r[i]));
        p->filename = "figs_data/finite_swarm3d_r"+std::to_string(r[i])+".txt";
        p->solve_fixed(250.,1e-6,true,500);
        
        printf("Finished with r=%g\n",r[i]);
        
    }
    
}

