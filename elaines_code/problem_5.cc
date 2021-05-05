//
//  problem_5.cpp
//  hw2_code
//
//  Created by Elaine Cunha on 3/3/21.
//

#include "problem_5.hh"
#include "problem_1_rk4d.hh"

/** The number of components in the ODE system. */
//const int ns=3;

int main() {
    
    double J=0.5, K=0.5, N=1250;
    
    p5_rk4d* o; o=new p5_rk4d(J, K, N);
//    o->solve_fixed(1,0.000001,false,0);
    o->solve_fixed(250.,0.000001,false,500);
    
}

