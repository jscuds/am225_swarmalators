#include <stdio.h>
#include "problem_1_rk4d.hh"
#include <iostream>
#include <cmath>
#include "custom_rng_code.h"
#include "swarm_ode.hh"

/** This class has functions to specify the ode for problem 5. */

// set up initial conditions
void swarm::init(double *q) {
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
void swarm::ff(double t_,double *in,double *out) {
    
    // for each agent
    for (int i=0;i<3*N;i+=3) {
        int idx_forcing = i/3;

        // initialize variables for agent i and start summations
        double x_i = in[i], y_i = in[i+1], theta_i = in[i+2];

        //forcing terms
        out[i] = vx[idx_forcing];
        out[i+1] = vy[idx_forcing];
        out[i+2] = w[idx_forcing];
        
        // for each pairwise interaction between agents
        for (int j=0;j<3*N;j+=3) {

            if (i != j) {
                bool eval_temp = eval_interaction(i,j,in);

                if (eval_temp == true){
                    // initialize variables for agent j
                    double x_j=in[j], y_j=in[j+1], theta_j=in[j+2];

                    // calculate intermediate terms
                    double dx = x_j - x_i;
                    double dy = y_j - y_i;
                    double xy_norm = sqrt(dx*dx + dy*dy);
                    double dtheta = theta_j - theta_i;
                    
                    // update output
                    out[i] += N_inv*((dx/xy_norm)*(1.+J*cos(dtheta))-(dx/(xy_norm*xy_norm)));
                    out[i+1] += N_inv*((dy/xy_norm)*(1.+J*cos(dtheta))-(dy/(xy_norm*xy_norm)));
                    out[i+2] += (K*N_inv)*sin(dtheta)/xy_norm;
                } 
            }
        }
    }
}

bool swarm::eval_interaction(int i,int j,double *in){
    if (r < 0){
        return true;
    } else if (checked_inter[i*N_int + j] == true){
        return interactions[i*N_int + j];
    } else{
        int idx1 = i*N_int + j;
        int idx2 = j*N_int + i;

        double x_j = in[j], y_j = in[j + 1];
        double x_i = in[i], y_i = in[j + 1];

        double r_temp = sqrt((x_i - x_j)*(x_i - x_j) + (y_i - y_j)*(y_i - y_j));

        if (r_temp < r){
            checked_inter[idx1] = true;
            checked_inter[idx2] = true;

            interactions[idx1] = true;
            interactions[idx2] = true;
            return true;
        } else{
            checked_inter[idx1] = true;
            checked_inter[idx2] = true;

            interactions[idx1] = false;
            interactions[idx2] = false;
            return false;
        }
    }
}