//
//  problem_5.hpp
//  hw2_code
//
//  Created by Elaine Cunha on 3/3/21.
//

#ifndef problem_5_hh
#define problem_5_hh

#include <stdio.h>
#include "problem_1_rk4d.hh"
#include <iostream>
#include <cmath>
#include "custom_rng_code.h"

/** This class has functions to specify the ode for problem 5. */
class p5_ode {
    public:
        // constructor
        const double J, K, N;
        double* const vx;
        double* const vy;
        double* const w;
        p5_ode(double J_, double K_, double N_): J(J_), K(K_), N(N_), vx(new double[int(N)]), vy(new double[int(N)]), w(new double[int(N)]) {
            // initialize w and v =0
            for (int i=0; i<N; i++) {
                vx[i]=0.01;
                vy[i]=0.01;
                w[i]=0.1;
            }
        }
        // destructor
        ~p5_ode() {
            delete [] vx;
            delete [] vy;
            delete [] w; 
        }
        // set up initial conditions
        void p5_ode_init(double *q) {
            double t_temp, r_temp, rand;
            custom_rng* c; c=new custom_rng(1);
            
            for (int i=0; i<3*N; i+=3) {
                // generate random number
                rand = c->doub();
                // get x and y positions
                t_temp = 2.*M_PI*rand;
                r_temp = c->doub();
                // initialize solution vector
                q[i] = r_temp*cos(t_temp);
                q[i+1] = r_temp*sin(t_temp);
                q[i+2] = 2.*M_PI*rand;
            }
            delete [] c;
        }
        // evaluate RHS of ODE system
        void p5_ode_ff(double t_,double *in,double *out) {
            
            // for each agent
            for (int i=0;i<3*N;i+=3) {
                // initialize variables for agent i and start summations
                double x_i=in[i], y_i=in[i+1], theta_i=in[i+2];
                out[i]=vx[i/3];
                out[i+1]=vy[i/3];
                out[i+2]=w[i/3];
                
                // for each pairwise interaction between agents
                for (int j=0;j<3*N;j+=3) {
                    if (i != j) {
                        // initialize variables for agent j
                        double x_j=in[j], y_j=in[j+1], theta_j=in[j+2];
                        
                        // calculate intermediate terms
                        double dx = x_j - x_i;
                        double dy = y_j - y_i;
                        double xy_norm = sqrt(dx*dx + dy*dy);
                        double dtheta = theta_j - theta_i;
                        
                        // update output
                        out[i] += (1./N)*((dx/xy_norm)*(1.+J*cos(dtheta))-(dx/(xy_norm*xy_norm)));
                        out[i+1] += (1./N)*((dy/xy_norm)*(1.+J*cos(dtheta))-(dy/(xy_norm*xy_norm)));
                        out[i+2] += (K/N)*sin(dtheta)/xy_norm;
                    }
                }
            }
        }
};

/** Class to solve the Brusselator problem with the adaptive fourth-order Runge-Kutta
 * method. */
class p5_rk4d : public p1_rk4d, public p5_ode {
    public:
        p5_rk4d(double J_, double K_, double N_) : p1_rk4d(3*N_), p5_ode(J_, K_, N_) {}
        virtual void ff(double t_,double *in,double *out) {
            p5_ode_ff(t_,in,out);
        }
        virtual void init() {p5_ode_init(q);}
};


#endif /* problem_5_hpp */
