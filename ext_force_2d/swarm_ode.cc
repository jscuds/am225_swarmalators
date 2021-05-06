#include <stdio.h>
#include <iostream>
#include <cmath>

#include "swarm_ode.hh"

/** This class has functions to specify the ode for problem 5. */

// set up initial conditions
void swarm::init(double *q) {
    double t_temp, r_temp, th0;
    custom_rng* c; c=new custom_rng(1);
    
    for (int i=0; i<3*N; i+=3) {

        // generate (r,th) for pt, sqrt-ing r
        t_temp = 2 * M_PI * (c->doub());
        r_temp = sqrt(c->doub());

        th0 = 2 * M_PI * (c->doub());

        // initialize solution vector
        q[i] = r_temp*cos(t_temp);
        q[i+1] = r_temp*sin(t_temp);
        q[i+2] = th0;
    }
    delete c;
}

// evaluate RHS of ODE system
void swarm::ff(double t_,double *in,double *out) {
    
    // for each agent
    for (int i=0;i<3*N;i+=3) {

        // initialize variables for agent i and start summations
        double x_i = in[i], y_i = in[i+1], theta_i = in[i+2];
        int idx_vw = i/3;

        // initial velocity terms
        if (F==0) { // without forcing
            out[i] = vx[idx_vw];
            out[i+1] = vy[idx_vw];
            out[i+2] = w[idx_vw];
        } else { // with forcing
            out[i] = 0.;
            out[i+1] = 0.;
            
            // add forcing to theta update
            double F_dx = F_locx - x_i;
            double F_dy = F_locy - y_i;
            double force_loc_norm = sqrt(F_dx*F_dx + F_dy*F_dy);
            out[i+2] = F*(cos(F_freq*t_ - theta_i)/force_loc_norm);
        }
        
        // for each pairwise interaction between agents
        for (int j=0;j<3*N;j+=3) {

            if (i != j) {
                bool eval_temp = eval_interaction(i/3,j/3,in);

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
    clear_checked();
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
            // std::cout << "i: " << i << "\n";
            // std::cout << "j: " << j << "\n";
            // std::cout << "idx1: " << idx1 << "\n";
            // std::cout << "idx2: " << idx2 << "\n";
            // std::cout << "size: " << N_int*N_int << "\n";
            checked_inter[idx1] = true;
            checked_inter[idx2] = true;

            interactions[idx1] = false;
            interactions[idx2] = false;
            return false;
        }
    }
}

void swarm::clear_checked(){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            checked_inter[i*N_int + j] = false;
        }
    }
}
