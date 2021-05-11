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
#include "point.hh"

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
        const int dim, N_int;
        double* const forcing;
        const double N_inv;
        bool* interactions;
        bool* checked_inter;
    
        // for grid of boxes
        const double n_boxes; // number of boxes in each dimension
        const int total_boxes, box_max; // total number of boxes, max number of points for each box
        double* num_boxes; // number of boxes along each dimension
        double* grid_dims; // tracks min/max dimensions for grid
        double* box_dims; // dxdydz dimensions of boxes
        int* pt_count; // array for storing total number of points in each box
        point** boxes; // pointer to array of boxes (each box is an array of points)


        swarm(double J_, double K_, int N_, int dim_, double r_): J(J_), K(K_), N(N_), r(r_), dim(dim_),
            N_int(int(N_)), forcing(new double[(dim + 1) * N_int]), N_inv(1./N),
            interactions(new bool[N_int * N_int]), checked_inter(new bool[N_int * N_int]),
            n_boxes(20.), total_boxes(int(pow(n_boxes,dim))), box_max(100),
            num_boxes(new double[dim]), grid_dims(new double[2*dim]), box_dims(new double[dim]),
            pt_count(new int[total_boxes]), boxes(new point *[total_boxes]) {
                
            // initialize w and v = 0
            for (int i=0; i<N; i++) {
                int idx = (dim + 1)*i;
                for (int j = 0; j < dim; j++){
                    forcing[idx + j] = 0.;
                }
                forcing[idx + dim] = 0.;

                for (int j = 0; j < N; j++){
                    checked_inter[i*N_int + j] = false;
                }
            }
                
            // initialize grid of boxes
            for (int i=0; i<dim; i++) num_boxes[i]=n_boxes;
            for (int i=0; i<total_boxes; i++) pt_count[i]=0;
            for (int i=0; i<total_boxes; i++) boxes[i]=new point[box_max];
                
        }
        // destructor
        ~swarm() {
            delete [] checked_inter;
            delete [] interactions;
            delete [] forcing;
            delete [] num_boxes;
            delete [] grid_dims;
            delete [] box_dims;
            delete [] pt_count;
            delete [] boxes;
        }
        // set up initial conditions
        void init(double *q);

        // evaluate RHS of ODE system
        void ff(double t_,double *in,double *out);
        void ff2(double t_,double *in,double *out);

        
    private:
    
        //get dimensions of grid
        void get_grid_dimensions(double *in, int dim);
    
        //populate boxes with points
        void populate_boxes(double *in);
    
        //clear boxes
        void clear_boxes();
    
        // evaluate equations for points i and j
        void evaluate_equations(int idx_i, int idx_j, double* x_i, double* x_j, double theta_i,
                                double theta_j, int dim, double *out, double *in);
        
        // extract xyz data from point
        void update_xyz(double *x, point p);
    
        // get high/low indices of boxes in subgrid
        void get_neighbors(double *x, double r, int *subgrid_idx, int dim);
    
        //eval finite distance cuttoff
        bool eval_interaction(int i,int j,double *in);

        //clear checked array (set all to false again)
        void clear_checked();

        //eval euclidean distance
        double euclidean_distance(double* x_i, double* x_j);
};

/** Class to solve the Brusselator problem with the adaptive fourth-order Runge-Kutta
 * method. */
class swarm_sol : public fsal_rk4d, public swarm {
    public:
        swarm_sol(double J_, double K_, double N_, double r_ = -1) : fsal_rk4d(3*N_), swarm(J_, K_, N_, 2, r_) {}
        void ff(double t_,double *in,double *out) override {
            swarm::ff(t_,in,out);
        }
        void init() override {
            swarm::init(q);
        }
};

class swarm_sol_3d : public fsal_rk4d, public swarm {
    public:
        swarm_sol_3d(double J_, double K_, double N_, double r_ = -1) : fsal_rk4d(4*N_), 
            swarm(J_, K_, N_, 3, r_) {}
        void ff(double t_,double *in,double *out) override {
            swarm::ff(t_,in,out);
        }
        void init() override {
            //todo:fix init for n dimensions
            swarm::init(q);
        }
};


#endif /* swarm_ode.hh */
