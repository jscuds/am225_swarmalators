//
//  problem_5.cpp
//  hw2_code
//
//  Created by Elaine Cunha on 3/3/21.
//
#include <cmath>
#include <memory>

#include "swarm_ode.hh"
#include "fsal_rk4d.hh"

/** The number of components in the ODE system. */
//const int ns=3;

int main() {
    
    double J = 0.1, K = 1., N = 250;
    double r = -1., F = 2., F_freq=3*M_PI/2, F_locx=0.1, F_locy=0.1;
    
//    std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,r));
//    o->filename = "unforce_test.txt";
//    o->solve_fixed(250,0.000001,true,500);
//
//    puts("Done with unforced test.");
    
    std::unique_ptr<swarm_sol> o2(new swarm_sol(J,K,N,r,F,F_freq,F_locx,F_locy));
    o2->filename = "force_test.txt";
    o2->solve_fixed(250,0.000001,true,500);
    
    puts("Done with forced test.");
    
//    std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,r,F,F_freq,F_locx,F_locy)); //,0.5));
//    o->filename = "FORCE_TEST_05_05_with_shift.txt";
//    o->solve_fixed(1,0.000001,true,5);

    //std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N)); //,0.5));
    //o->filename = "ORIG_05_05.txt";
    //o->solve_fixed(250,0.000001,false,5);


    // Helen's finite test
    //std::unique_ptr<swarm_sol> o(new swarm_sol(J,K,N,0.5));
    //swarm_sol* o; o=new swarm_sol(J,K,N,0.5);
    //o->filename = "FINITE_TEST.txt";
    //o->solve_fixed(10.,0.000001,true,5);

    //delete o;
    
}

