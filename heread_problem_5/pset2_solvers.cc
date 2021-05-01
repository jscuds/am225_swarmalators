#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <utility>
#include <cstdio>
#include <iomanip>
#include "pset2_odes.hh"

typedef std::pair<int, double> sol_out;
typedef std::pair<int, double*> count_fin;

//euclidean distance
double euclidean_distance(double *y1, double *y2, int n){
    double summed_temp = 0;
    for (int i = 0; i < n; i++){
        summed_temp += (y1[i] - y2[i]) * (y1[i] - y2[i]);
    }
    return std::sqrt(summed_temp);
}

//fsal solver
double* rk4_final_check(double tmax, int steps, ode *problem){
    double h = tmax/steps;
    int n = problem->n;

    double *y = new double[n];
    double t = 0;

    double *k1 = new double[n];
    double *k2 = new double[n];
    double *k3 = new double[n];
    double *k4 = new double[n];
    double *dq = new double[n];

    problem->initialize(y);

    for (int i = 0; i < steps; i++){
        //k1
        problem->equation(t,y,k1);

        //k2
        for (int i = 0; i < n; i++) dq[i] = y[i] + h*k1[i]/2;
        problem->equation(t + h/2,dq,k2);

        //k3
        for (int i = 0; i < n; i++) dq[i] = y[i] + h*k2[i]/2;
        problem->equation(t + h/2,dq,k3);

        //k4
        for (int i = 0; i < n; i++) dq[i] = y[i] + h*k3[i];
        problem->equation(t + h,dq,k4);

        for (int i = 0; i < n; i++){
            y[i] = y[i] + h/6 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        }
        //update
        t += h;

    }
    delete [] dq;
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    
    return y;
}

count_fin fsal_solver(double tmax, double lambda, ode *problem, bool writeq, std::string filename = ""){

    //set values for step size selection
    double facmin = 1/3;
    double facmax = 3;
    double fac = 0.9;

    //initial h
    double h = 0.01;

    //dim of problem
    int n = problem->n;
    double m = n;
    //std::cout << n << "\n";

    //create y and dy & initialize t
    double *y = new double[n];
    double *y1 = new double[n];
    double *y1_hat = new double[n];
    double *sc = new double[n];
    double t = 0;

    problem->initialize(y);
    //todo: get rid of this
    std::cout << "initialized" << "\n";

    //actual solver code
    double *k1 = new double[n];
    double *k2 = new double[n];
    double *k3 = new double[n];
    double *k4 = new double[n];
    double *k5 = new double[n];
    double *dq = new double[n];

    double error; //make error
    double error_temp_sum = 0;

    //make output file, only write to it if bool writeq is true (saves time if false)
    std::ofstream matrix_file(filename, std::ofstream::trunc);
    if (writeq){
        matrix_file << 0 << ",";
    }
    
    for (int i = 0; i < n; i++){
        if (writeq){
            matrix_file << y[i];
            if (i == n-1){
                matrix_file << std::endl;
            } else {
                matrix_file << ",";
            }
        }
    }

    problem->equation(0,y,k1);
    int func_calls = 1;
    bool flag_last_call = false;

    while (t < tmax || flag_last_call == true){
        //time (todo: delete this)
        //std::cout <<"t: " << t << "\n";
        //k2
        for (int i = 0; i < n; i++) dq[i] = y[i] + (h/3 * k1[i]);
        problem->equation(t + h/3, dq, k2);
        func_calls += 1;

        //k3
        for (int i = 0; i < n; i++) dq[i] = y[i] + h * (-1./3 * k1[i] + k2[i]);
        problem->equation(t + 2*h/3, dq, k3);
        func_calls += 1;

        //k4
        for (int i = 0; i < n; i++) dq[i] = y[i] + h * (k1[i] - k2[i] + k3[i]);
        problem->equation(t + h, dq, k4);
        func_calls += 1;

        //k5
        for (int i = 0; i < n; i++) dq[i] = y[i] + h/8 * (k1[i] + 3*k2[i] + 3*k3[i] + k4[i]);
        problem->equation(t + h, dq, k5);
        func_calls += 1;
        
        for (int i = 0; i < n; i++){
            y1[i] = dq[i];
            y1_hat[i] = y[i] + h/12 * (k1[i] + 6*k2[i] + 3*k3[i] + 2*k5[i]);
            sc[i] = lambda + lambda * std::max(std::abs(y[i]), std::abs(y1[i]));
            error_temp_sum += ((y1[i] - y1_hat[i])/sc[i]) * ((y1[i] - y1_hat[i])/sc[i]);
            //std::cout << "y1: " << y1[i] <<", "<<"y1hat: " << y1_hat[i]<<"\n"; 
        }

        error = std::sqrt(1/m * error_temp_sum);
        //std::cout<<"error: " << error<<"\n";
        //should reject step if error is too large
        if (error < 1){
            //writing:
            if (writeq){
                matrix_file << t + h << ",";
                for (int i = 0; i < n; i++){
                    matrix_file << y1[i];
                    if (i == n-1){
                        matrix_file << std::endl;
                    } else {
                        matrix_file << ",";
                    }
                    
                }
            }
            //setting y and k1 for next round
            for (int i = 0; i < n; i++){
                y[i] = y1[i];
                k1[i] = k5[i];
            }
            t += h;
        } else {
            flag_last_call = false;
        }
        h = h * std::min(facmax, std::max(facmin, fac * pow(error, -0.25)));
        //std::cout <<"h: " << h << "\n";
        error_temp_sum = 0;

        if (t + h > tmax && flag_last_call == false) {
            h = tmax - t; //should be less than old h
            flag_last_call = true; //should be last call
        } else if (flag_last_call) {
            flag_last_call = false;
        }
    }
    /*if (writeq){
        matrix_file.close();
    }*/
    matrix_file.close();

    // Free dynamically allocated memory
    delete [] dq;
    delete [] k5;
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] sc;
    delete [] y1_hat;
    delete [] y1;

    count_fin fin_return(func_calls, y);
    return fin_return;
}

/*if (flag_last_call == true) {
                //reference soln below is for bruss for t_max = 20
                //double reference_soln[2] = {0.4986370712683506, 4.596780349452039};
                //todo: fix this for p5
                //double* reference_soln = rk4_final_check(tmax, 100000, problem);
                //final_error = euclidean_distance(y1,reference_soln,n);
                final_error = 1.0;

                //delete [] reference_soln;
}*/

//prob 1 processing
sol_out prob1_processing(double tmax, double lambda, ode *problem, bool writeq, std::string filename = ""){
    count_fin final_result = fsal_solver(tmax,lambda,problem,writeq,filename);
    int n = problem->n;
    double* reference_soln = rk4_final_check(tmax, 100000, problem);
    double final_error = euclidean_distance(final_result.second,reference_soln,n);
    sol_out prob1_res(final_result.first,final_error);
    return prob1_res;
}

//fsal solver final only
count_fin fsal_prob7(double tmax, double atol, double rtol, ode *problem){

    //set values for step size selection
    double facmin = 1/3;
    double facmax = 3;
    double fac = 0.9;

    //initial h
    double h = 0.01;

    //dim of problem
    int n = problem->n;
    double m = n;

    //create y and dy & initialize t
    double *y = new double[n];
    double *y1 = new double[n];
    double *y1_hat = new double[n];
    double *sc = new double[n];
    double t = 0;

    problem->initialize(y);

    //actual solver code
    double *k1 = new double[n];
    double *k2 = new double[n];
    double *k3 = new double[n];
    double *k4 = new double[n];
    double *k5 = new double[n];
    double *dq = new double[n];

    double error; //make error
    double error_temp_sum = 0;

    problem->equation(0,y,k1);
    int func_calls = 1;
    bool flag_last_call = false;
    double final_error;

    while (t < tmax || flag_last_call == true){
        //k2
        for (int i = 0; i < n; i++) dq[i] = y[i] + (h/3 * k1[i]);
        problem->equation(t + h/3, dq, k2);
        func_calls += 1;

        //k3
        for (int i = 0; i < n; i++) dq[i] = y[i] + h * (-1./3 * k1[i] + k2[i]);
        problem->equation(t + 2*h/3, dq, k3);
        func_calls += 1;

        //k4
        for (int i = 0; i < n; i++) dq[i] = y[i] + h * (k1[i] - k2[i] + k3[i]);
        problem->equation(t + h, dq, k4);
        func_calls += 1;

        //k5
        for (int i = 0; i < n; i++) dq[i] = y[i] + h/8 * (k1[i] + 3*k2[i] + 3*k3[i] + k4[i]);
        problem->equation(t + h, dq, k5);
        func_calls += 1;
        
        for (int i = 0; i < n; i++){
            y1[i] = dq[i];
            y1_hat[i] = y[i] + h/12 * (k1[i] + 6*k2[i] + 3*k3[i] + 2*k5[i]);
            sc[i] = atol + rtol * std::max(std::abs(y[i]), std::abs(y1[i]));
            error_temp_sum += ((y1[i] - y1_hat[i])/sc[i]) * ((y1[i] - y1_hat[i])/sc[i]);
        }

        error = std::sqrt(1/m * error_temp_sum);
        //should reject step if error is too large
        if (error < 1){
            //setting y and k1 for next round
            for (int i = 0; i < n; i++){
                y[i] = y1[i];
                k1[i] = k5[i];
            }
            t += h;
        } else {
            flag_last_call = false;
        }
        h = h * std::min(facmax, std::max(facmin, fac * pow(error, -0.25)));
        error_temp_sum = 0;

        if (t + h > tmax && flag_last_call == false) {
            h = tmax - t; //should be less than old h
            flag_last_call = true; //should be last call
        } else if (flag_last_call) {
            flag_last_call = false;
        }
    }

    // Free dynamically allocated memory
    delete [] dq;
    delete [] k5;
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] sc;
    delete [] y1_hat;
    delete [] y1;

    count_fin fin_return(func_calls, y);
    return fin_return;
}