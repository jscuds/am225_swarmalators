//
//  problem_5.hpp
//  hw2_code
//
//  Created by Elaine Cunha on 3/3/21.
//

#ifndef swarm_ode_hh
#define swarm_ode_hh

#include <stdio.h>
#include <iostream>
#include <cmath>

#include "fsal_rk4d.hh"

/** Custom random number generator based of "Ran" routine in Numerical Recipes
 * by Press et al. */
class custom_rng {
    public:
        unsigned long a,b,c;
        custom_rng(unsigned long seed) : b(4101842887655102017L), c(1) {
            if(sizeof(unsigned long)<8) {
                fputs("Error: 'unsigned long' type too short\n",stderr);
                exit(1);
            }
            a=seed^b;int64();
            b=a;int64();
            c=b;int64();
        }
        unsigned long int64() {
            a=a*2862933555777941757L+7046029254386353087L;
            b^=b>>17;b^=b<<31;b^=b>>8;
            c=4294957665U*(c&0xffffffff)+(c>>32);
            unsigned long d=a^(a<<21);
            d^=d>>35;d^=d<<4;
            return (d+b)^c;
        }
        inline double doub() {
            return 5.42101086242752217E-20*int64();
        }
};

/** This class has functions to specify the ode for problem 5. */
class swarm {
    public:
        // constructor
        const double J, K, N, r;
        double* const vx;
        double* const vy;
        double* const w;
        const double N_inv;
        const int N_int;
        bool* interactions;
        bool* checked_inter;


        swarm(double J_, double K_, double N_, double r_): J(J_), K(K_), N(N_), r(r_),
            vx(new double[int(N)]), vy(new double[int(N)]), w(new double[int(N)]), N_inv(1./N),
            N_int(int(N)), interactions(new bool[N_int * N_int]), checked_inter(new bool[N_int * N_int]) {
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
            delete [] checked_inter;
            delete [] interactions;
            delete [] w;
            delete [] vy;
            delete [] vx;
        }
        // set up initial conditions
        void init(double *q);

        // evaluate RHS of ODE system
        void ff(double t_,double *in,double *out);

        
    private:
        //eval finite distance cuttoff
        bool eval_interaction(int i,int j,double *in);

        //clear checked array (set all to false again)
        void clear_checked();
};

/** Class to solve the Brusselator problem with the adaptive fourth-order Runge-Kutta
 * method. */
class swarm_sol : public fsal_rk4d, public swarm {
    public:
        swarm_sol(double J_, double K_, double N_, double r_ = -1) : fsal_rk4d(3*N_), swarm(J_, K_, N_, r_) {}
        void ff(double t_,double *in,double *out) override {
            swarm::ff(t_,in,out);
        }
        void init() override {
            swarm::init(q);
        }
};


#endif /* swarm_ode.hh */
