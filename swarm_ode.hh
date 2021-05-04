//
//  problem_5.hpp
//  hw2_code
//
//  Created by Elaine Cunha on 3/3/21.
//

#ifndef swarm_ode_hh
#define swarm_ode_hh

#include <stdio.h>
#include "problem_1_rk4d.hh"
#include <iostream>
#include <cmath>
#include "custom_rng_code.h"

/** This class has functions to specify the ode for problem 5. */
class swarm {
    public:
        // constructor
        const double J, K, N, N_inv, r;
        double* const vx;
        double* const vy;
        double* const w;
        const int N_int;
        bool* interactions;
        bool* checked_inter;


        swarm(double J_, double K_, double N_, double r): J(J_), K(K_), N(N_), r(r_),
            vx(new double[int(N)]), vy(new double[int(N)]), w(new double[int(N)]), N_inv = 1./N,
            interactions(new bool[int(N)*int(N)]), N_int(int(N)), checked_inter(new bool[N_int*N_int]) {
            // initialize w and v = 0
            for (int i=0; i<N; i++) {
                vx[i]=0.01;
                vy[i]=0.01;
                w[i]=0.1;
                for (int j = 0; j < N; j++){
                    checked_inter[i*N_int + j] = false;
                }
            }
        }
        // destructor
        ~swarm() {
            delete [] vx;
            delete [] vy;
            delete [] w; 
        }
        // set up initial conditions
        void init(double *q);

        // evaluate RHS of ODE system
        void ff(double t_,double *in,double *out);

        
    private:
        //eval finite distance cuttoff
        bool eval_interaction(int i,int j,double *in);
};

/** Class to solve the Brusselator problem with the adaptive fourth-order Runge-Kutta
 * method. */
class p5_rk4d : public p1_rk4d, public swarm {
    public:
        p5_rk4d(double J_, double K_, double N_, double r_ = -1) : p1_rk4d(3*N_), swarm(J_, K_, N_, r_) {}
        virtual void ff(double t_,double *in,double *out) {
            ff(t_,in,out);
        }
        virtual void init() {init(q);}
};


#endif /* problem_5_hpp */
