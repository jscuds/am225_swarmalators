#include <stdio.h>
#include <iostream>
#include <cmath>
#include <stdexcept>

#include "swarm_ode.hh"
#include "point.hh"

/** This class has functions to specify the ode. */

// set up initial conditions
void swarm::init(double *q) {
    if (dim == 2){
        init2d(q);
    } else{
        custom_rng* c; c = new custom_rng(1);
        double* state_test = new double[dim];

        //make array of all zeros to pass in later
        double* zero_vec = new double[dim];
        for (int k = 0; k < dim; k++) zero_vec[k] = 0;

        for (int i = 0; i < N_int; i++){
            int idx_i = (dim + 1) * i;
            double L = 2;
            while (L > 1){
                for (int j = 0; j < dim; j++){
                    state_test[j] = 2*(c->doub()) - 1;
                }
                L = euclidean_distance(state_test,zero_vec);
            }
            //if it passed, update q
            for (int j = 0; j < dim; j++){
                q[idx_i + j] = state_test[j];
            }
            q[idx_i + dim] = 2*M_PI*(c->doub());
        }
        delete [] state_test;
        delete [] zero_vec;
        delete c;
    }
//    std::cout << "initialized" << "\n";
}

void swarm::init2d(double* q){
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

// determine min and max dimensions for x,y,z
void swarm::get_grid_dimensions(double *in, int dim) {
    
    // initialize grid_dims
    for (int i = 0;i < dim;i++) {
        grid_dims[i*2] = in[i]; //min
        grid_dims[i*2+1] = in[i]; //max
    }
    
    // find min and max dims across all points
    int idx_i;
    for (int i = 0;i < N;i++) {
        idx_i = (dim + 1)*i;
        for (int j = 0;j < dim;j++) { //for x,y,z
            if (grid_dims[j*2] > in[idx_i+j]) grid_dims[j*2] = in[idx_i+j]; //update min
            if (grid_dims[j*2+1] < in[idx_i+j]) grid_dims[j*2+1] = in[idx_i+j]; //update max
        }
    }

    // update box dimensions (dx, dy, dz for each box)
    for (int i = 0;i < dim;i++) {
        box_dims[i] = fabs(grid_dims[i+1] - grid_dims[i])/num_boxes[i];
    }
}

// put points into boxes
void swarm::populate_boxes(double *in) {
    
    // initialize variables
    int equiv_box_idx[dim];
    int idx_i, box_idx, count;
    
    // iterate over every point
    for (int i = 0;i < N;i++) {
        idx_i = (dim + 1)*i;
        
        // convert from (x,y) to (i,j) of box in grid
        for (int j = 0;j < dim;j++) { // shift grid so in positive domain before converting to i,j
            equiv_box_idx[j] = int((in[idx_i+j] - grid_dims[j*2])/box_dims[j]);
            // account for overflow (sometimes eq. above results in idx = num_boxes, so move to
            // closest box inside grid)
            if (equiv_box_idx[j] >= num_boxes[j]) equiv_box_idx[j] = num_boxes[j]-1;
        }
                
        // convert (i,j) of box in grid to idx within boxes array
        box_idx = equiv_box_idx[0]*num_boxes[0] + equiv_box_idx[1];
        if (dim==3) box_idx += equiv_box_idx[2]*(num_boxes[0]*num_boxes[1]);
        
        // add to count and check if max points exceeded for box
        count = pt_count[box_idx];
        if (count+1 > box_max) {
            throw std::invalid_argument( "Max points per box exceeded." );
        }
        
        // add point to correct box array
        if (dim==2) {
            boxes[box_idx][count].set(in[idx_i], in[idx_i+1], 0., in[idx_i+2], idx_i);
            pt_count[box_idx] += 1;
        } else {
            boxes[box_idx][count].set(in[idx_i], in[idx_i+1], in[idx_i+2], in[idx_i+3], idx_i);
            pt_count[box_idx] += 1;
        }
    }
}

// clear boxes, reset point count
void swarm::clear_boxes() {
    
    for (int i = 0;i < total_boxes;i++) {
        for (int j = 0;j < pt_count[i];j++) boxes[i][j].clear(); //clear data
        pt_count[i] = 0; // reset point count for box
    }
}

// extract xyz data from point
void swarm::update_xyz(double *x, point p) {
    x[0] = p.x;
    x[1] = p.y;
    x[2] = p.z;
}

// get high/low indices of boxes in subgrid
void swarm::get_neighbors(double* x, double r, int *subgrid_idx, int dim) {
    
    // REFERENCE:
    // num_boxes = # boxes along each dimension
    // box_dims = [dx, dy, dz] of boxes
    // grid_dims = min/max dims of grid [xmin, xmax, ymin, ymax, ...]
    // subgrid_idx = [i_min, i_max, j_min, j_max, ...]
    
    double low, high;
    int ilow, ihigh;
    
    if (r<0) { // if infinite radius, set subgrid to include all boxes
        for (int i = 0;i < dim;i++) {
            subgrid_idx[i*2] = 0;
            subgrid_idx[i*2+1] = num_boxes[i]-1;
        }
    } else { // calculate min/max box indices
        for (int i = 0;i < dim;i++) {
            // get min and max coords of subgrid for dimension i
            low = fmax(grid_dims[i*2], x[i] - r);
            high = fmin(grid_dims[i*2+1], x[i] + r);
        
            // convert to box indices in i,j,k coords
            ilow = int((low - grid_dims[i*2])/box_dims[i]); //shift grid so in positive domain
            ihigh = int((high - grid_dims[i*2])/box_dims[i]); //shift grid so in positive domain
            
            // account for overflow (sometimes eq. above results in idx = num_boxes, so move to
            // closest box inside grid)
            if (ihigh >= num_boxes[i]) ihigh = num_boxes[i]-1;
        
            // update subgrid_idx
            subgrid_idx[i*2] = ilow;
            subgrid_idx[i*2+1] = ihigh;
        }
    }
}

void swarm::evaluate_equations(int idx_i, int idx_j, double* x_i, double* x_j, double theta_i,
                               double theta_j, int dim, double *out, double *in) {
    
    double dx;
    bool eval_temp = eval_interaction(idx_i,idx_j,in);
    
    if (eval_temp == true){

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

// evaluate RHS of ODE system
void swarm::ff(double t_,double *in,double *out) {
    
    double* x_i; x_i = new double[3];
    double* x_j; x_j = new double[3];
    double theta_i,theta_j;
    int* subgrid_idx; subgrid_idx = new int[2*dim];
    int idx_i, idx_j, box_idx;
    point p_i, p_j;
    
    // initialize grid with boxes
    get_grid_dimensions(in,dim);
    clear_boxes();
    
    // fill in boxes with data from in array
    populate_boxes(in);
    
    // iterate through points and evaluate equations
    // for box in grid
    for (int i = 0;i < total_boxes;i++) {
        // for point in box
        for (int p = 0;p < pt_count[i];p++) {
            // unpack point elements
            p_i = boxes[i][p];
            update_xyz(x_i, p_i);
            theta_i = p_i.theta;
            idx_i = p_i.idx;
            
            // initialize forcing
            for (int j = 0; j < dim + 1;j++){
                out[idx_i+j] = forcing[idx_i+j];
            }
            
            // calculate subgrid
            get_neighbors(x_i, r, subgrid_idx, dim);
            
            // iterate through subgrid boxes
            for (int i_box = subgrid_idx[0];i_box <= subgrid_idx[1];i_box++) { //width of subgrid
                for (int j_box = subgrid_idx[2];j_box <= subgrid_idx[3];j_box++) { //height of subgrid
                    if (dim==2) {
                        box_idx = i_box*num_boxes[0] + j_box;
                        // iterate through points in subgrid boxes
                        for (int pp = 0;pp < pt_count[box_idx];pp++) {
                            // unpack point elements
                            p_j = boxes[box_idx][pp];
                            update_xyz(x_j, p_j);
                            theta_j = p_j.theta;
                            idx_j = p_j.idx;
                            
                            // update equations
                            if (idx_i!=idx_j) {
                                evaluate_equations(idx_i, idx_j, x_i, x_j, theta_i, theta_j, dim, out, in);
                            }
                        }
                    } else {
                        for (int k_box = subgrid_idx[4];k_box <= subgrid_idx[5];k_box++) {
                            box_idx = i_box*num_boxes[0] + j_box + k_box*(num_boxes[0]*num_boxes[1]);
                            // iterate through points in subgrid boxes
                            for (int pp = 0;pp < pt_count[box_idx];pp++) {
                                // unpack point elements
                                p_j = boxes[box_idx][pp];
                                update_xyz(x_j, p_j);
                                theta_j = p_j.theta;
                                idx_j = p_j.idx;
                                
                                // update equations
                                if (idx_i!=idx_j) {
                                    evaluate_equations(idx_i, idx_j, x_i, x_j, theta_i, theta_j, dim, out, in);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    clear_checked();
    delete [] x_i;
    delete [] x_j;
    delete [] subgrid_idx;
    
//    for (int z=0;z<3*N;z++) printf("%g\n",in[z]);
//    if (true) {
//        throw std::invalid_argument( "STOP." );
//    }
}


// evaluate RHS of ODE system
void swarm::ff2(double t_,double *in,double *out) {
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
                
//                bool eval_temp = eval_interaction(i,j,in);
                idx_j = (dim + 1) * j;
                bool eval_temp = eval_interaction(idx_i,idx_j,in);
                
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
            idx_j = (dim + 1) * j;
        }
    }
//    clear_checked();
//    for (int z=0;z<3*N;z++) printf("%g\n",out[z]);
//    if (true) {
//     throw std::invalid_argument( "STOP." );
//    }
}


void swarm::clear_checked(){
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            checked_inter[i*N_int + j] = false;
        }
    }
}

bool swarm::eval_interaction(int idx_i,int idx_j,double *in){
    
    // convert idx_i and idx_j to indices of *in
    int i = idx_i/(dim+1);
    int j = idx_j/(dim+1);
    
    if (r < 0){
        return true;
    } else if (checked_inter[i*N_int + j] == true){
        return interactions[i*N_int + j];
    } else{
        int idx1 = i*N_int + j;
        int idx2 = j*N_int + i;

//        int idx_i = (dim + 1)*i;
//        int idx_j = (dim + 1)*j;
        

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


