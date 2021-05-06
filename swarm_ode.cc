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
    //std::cout << "t: " << t_ << "\n";
    
    double* x_i;
    double* x_j;
    double theta_i,theta_j,dx;
    int idx_i, idx_j;

    // for each agent
    for (int i = 0;i < N;i++) {
        idx_i = (dim + 1)*i;
        
        // initialize variables for agent i and start summations
        // xx_i = in[idx_i] = *(in + idx_i)
        x_i = in + idx_i;
        theta_i = in[idx_i + dim];

        //forcing terms--[note/todo: should we just initially copy the whole forcing to out?]
        for (int j = 0; j < dim + 1; j++){
            out[idx_i + j] = forcing[idx_i + j];

        }
        
        // for each pairwise interaction between agents
        for (int j = 0;j < N;j++) {

            if (i != j) {
                bool eval_temp = eval_interaction(i,j,in);

                if (eval_temp == true){
                    // initialize variables for agent j
                    idx_j = (dim + 1) * j;
                    x_j = in + idx_j;
                    theta_j = in[idx_j + dim];

                    // calculate intermediate terms
                    double L = euclidean_distance(x_i,x_j);
                    double L_mult = L;
                    double dtheta = theta_j - theta_i;
                    
                    for (int k = 1; k < dim; k++){
                        L_mult *= L;
                        //this is: L_mult = L_mult * L;
                        //should return L^2 for dim = 2, L^3 for dim = 3, etc
                    }
                    double L_inv = 1.0/L;
                    double L_fac = 1.0/L_mult;

                    // update output
                    for (int k = 0; k < dim; k++){
                        dx = x_j[k] - x_i[k];
                        out[idx_i + k] += N_inv * (dx*L_inv*(1 + J*cos(dtheta)) - dx*L_fac);
                    }
                    out[idx_i + dim] += (K*N_inv)*sin(dtheta) * L_inv;
                } 
            }
        }
    }
    clear_checked();
}

void swarm::clear_checked(){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            checked_inter[i*N_int + j] = false;
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

        int idx_i = (dim + 1)*i;
        int idx_j = (dim + 1)*j;

        double r_temp = euclidean_distance(in+idx_i,in+idx_j);

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

double swarm::euclidean_distance(double* x_i, double* x_j){
    double summed_temp = 0;
    for (int idx = 0; idx < dim; idx++){
        summed_temp += (x_i[idx] - x_j[idx]) * (x_i[idx] - x_j[idx]);
    }
    return sqrt(summed_temp);
}


