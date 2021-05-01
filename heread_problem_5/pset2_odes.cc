#include <cmath>
#include <cstdlib>
#include "pset2_odes.hh"
#include "pset2_solvers.hh"

//bruss
void brusselator::initialize(double *y){
    y[0] = 1.5;
    y[1] = 3;
}
void brusselator::equation(double t, double *y, double *dy){
    dy[0] = 1 + (y[0] * y[0] * y[1]) - 4*y[0];
    dy[1] = 3*y[0] - (y[0] * y[0] * y[1]);
}

//osc_1b
void osc_1b::initialize(double *y){
    y[0] = 1;
    y[1] = 0;
}
void osc_1b::equation(double t, double *y, double *dy){
    dy[0] = -t * y[1];
    dy[1] = t * y[0];
}

//dis_7
void dis_7::initialize(double *y){
    y[0] = 1;
    y[1] = 0;
}
void dis_7::equation(double t, double *y, double *dy){
    if (std::abs(y[0]) < std::abs(y[1])){
        dy[0] = -y[1];
        dy[1] = 0;
    } else {
        dy[0] = 0;
        dy[1] = y[0];
    }
}

//agents
inline double rnd(){
    return (1.0/RAND_MAX) * ((double) rand());
}
//rnd gives number between (0,1)
void agents::initialize(double *y){
    srand(0);
    int idx_start;
    for (int i = 0; i < N; i++){
        idx_start = 3*i;
        //unit disk (rip)
        double polar_r = sqrt(rnd());
        double polar_th = 2*M_PI*rnd();
        double x0 = polar_r * cos(polar_th);
        double y0 = polar_r * sin(polar_th);
        double th0 = 2*M_PI*rnd();
        double initial[] = {x0,y0,th0};
        for (int j = 0; j < 3; j++) y[idx_start+j] = initial[j];
        //M_PI
    }
}
void agents::equation(double t, double *y, double *dy){
    for (int i = 0; i < N; i++){
        int idx_start = i * 3;
        double temp_sum[] = {0,0};
        double temp_sum_th = 0;
        double* x_i = y + idx_start;
        double* theta_i = y + idx_start + 2;
        double* x_j; 
        double* theta_j;
        double length;
        int j_idx_start;
        for (int j = 0; j < N; j++){
            if (j != i){
                j_idx_start = j * 3;
                x_j = y + j_idx_start;
                theta_j = y + j_idx_start + 2;
                length = euclidean_distance(x_i,x_j,2);
                for (int k = 0; k < 2; k++){
                    temp_sum[k] += (x_j[k] - x_i[k])/length * (1 + J*cos(theta_j[0] - theta_i[0]));
                    temp_sum[k] -= (x_j[k] - x_i[k])/(length * length);
                }
                temp_sum_th += sin(theta_j[0] - theta_i[0])/length;
            }

        }
        for (int j = 0; j < 2; j++) dy[idx_start + j] = 1./N * temp_sum[j];
        dy[idx_start + 2] = K/N * temp_sum_th;
    }
}

//agents 3d
void agents3d::initialize(double *y){
    srand(0);
    int idx_start;
    for (int i = 0; i < N; i++){
        idx_start = 4*i;
        double length = 3;
        double x0,y0,z0;
        while (length > 1){
            x0 = 2*rnd() - 1;
            y0 = 2 * rnd() - 1;
            z0 = 2 * rnd() - 1;
            length = sqrt(x0*x0 + y0*y0 + z0*z0);
        }
        
        //unit sphere (rip)
        // double polar_r = pow(rnd(),1/3);
        // double polar_th = 2*M_PI*rnd();
        // double polar_phi = M_PI * rnd();
        // double x0 = polar_r * sin(polar_phi) * cos(polar_th);
        // double y0 = polar_r * sin(polar_phi) * sin(polar_th);
        // double z0 = polar_r * cos(polar_phi);
        double th0 = 2*M_PI*rnd();
        double initial[] = {x0,y0,z0,th0};
        for (int j = 0; j < 4; j++) y[idx_start+j] = initial[j];
        //M_PI
    }
}
void agents3d::equation(double t, double *y, double *dy){
    for (int i = 0; i < N; i++){
        int idx_start = i * 4;
        double temp_sum[] = {0,0,0};
        double temp_sum_th = 0;
        double* x_i = y + idx_start;
        double* theta_i = y + idx_start + 3;
        double* x_j; 
        double* theta_j;
        double length;
        int j_idx_start;
        for (int j = 0; j < N; j++){
            if (j != i){
                j_idx_start = j * 4;
                x_j = y + j_idx_start;
                theta_j = y + j_idx_start + 3;
                length = euclidean_distance(x_i,x_j,3);
                for (int k = 0; k < 3; k++){
                    temp_sum[k] += (x_j[k] - x_i[k])/length * (1 + J*cos(theta_j[0] - theta_i[0]));
                    temp_sum[k] -= (x_j[k] - x_i[k])/(length * length * length);
                }
                temp_sum_th += sin(theta_j[0] - theta_i[0])/length;
            }

        }
        for (int j = 0; j < 3; j++) dy[idx_start + j] = 1./N * temp_sum[j];
        dy[idx_start + 3] = K/N * temp_sum_th;
    }
}