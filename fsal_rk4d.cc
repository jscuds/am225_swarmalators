#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>

#include "fsal_rk4d.hh"

/** Initializes the fourth-order Runge-Kutta solver, allocating memory and
 * setting constants.
 * \param[in] dof_ the number of degrees of freedom. */
fsal_rk4d::fsal_rk4d(int dof_) : dof(dof_), fcount(0), t(0.), q(new double[dof]),
    dq(new double[dof]), dq_hat(new double[dof]),
    k1(new double[dof]), k2(new double[dof]), k3(new double[dof]),
    k4(new double[dof]), k5(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
fsal_rk4d::~fsal_rk4d() {
    delete [] k5;
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq_hat;
    delete [] dq;
    delete [] q; 
}

/** Prints the current state of the solution. */
void fsal_rk4d::print(double t_, double *in) {
    //todo: change to output to fxn

    //std::ofstream file_out(filename, std::ofstream::trunc);
    std::ofstream file_out(filename,std::ofstream::app);

    file_out << t_ << " ";
    for(int i = 0;i < dof; i++) {
        file_out << in[i] << " ";
    }
    file_out << "\n";

    file_out.close();
}

void fsal_rk4d::print_state(double t_, double *in) {
    //like print but doesn't print time
    std::ofstream file_out(filename,std::ofstream::app);

    for(int i = 0;i < dof; i++) {
        file_out << in[i] << " ";
    }
    file_out << "\n";
    file_out.close();
}

/** Solves the ODE problem with a fixed integration step.
 * \param[in] duration the time to solve for.
 * \param[in] steps the number of integration steps
 * \param[in] output whether to print each integration step. */
void fsal_rk4d::solve_fixed(double duration,double Atol,bool output,int d_steps) {

    // Set up initial condition and compute timestep. Use the t_den variable to
    // mark where the dense output has got to.
    init();

    if (output){
        std::ofstream file_out(filename, std::ofstream::trunc);
        file_out.close();
    }

    double dt = 0.01,t_den = 0.,dt_den;
    double dt_new = 0.01;
    bool last = false;
    int num_steps = 0;

    // Perform integration steps
    if(output) print(t,q);
    if(d_steps>0) {
        dt_den=duration/d_steps;
        print(t,q);
    }
   
    ff(t,q,k1);
    while(not last) {
        if (t+dt > duration) {
            dt = duration - t;
            last=true;
        }
        
        step(dt);
        
        // Calculate error
        double err=0., err_temp=0., sc_i=0.;
        for(int i=0;i<dof;i++) {
            sc_i = Atol + Atol*std::max(abs(q[i]),abs(dq[i]));
            err_temp += pow((dq[i]-dq_hat[i])/sc_i,2.);
        }
        err = sqrt(err_temp/dof);

        // Calculate h_new
        dt_new = dt*std::min(3., std::max(1./3., 0.9*pow(1./err, 0.25)));
        
        // Decide whether to take step
        if (err>1. and last) {// last step but step is rejected
            dt=dt_new;
            last=false;
        } else if (err>1. and not last) {// rejected, but not last step
            dt=dt_new;
        } else { // take step
            t += dt;
            num_steps += 1;

            // Do any dense output interpolation
            if(d_steps>0) {
                while(t_den+dt_den<t) {
                    t_den += dt_den;
                    dense_output(1.+(t_den-t)/dt,dt);
                    if ((fabs(t_den-10.) < 1e-12) or (fabs(t_den-20.) < 1e-12) or 
                        (fabs(t_den-50.) < 1e-12) or (fabs(t_den-200.) < 1e-12)) {
                        print(t_den,k3);
                    }
                }
            }
            
            dt = dt_new;
            memcpy(q,dq,dof*sizeof(double));
            memcpy(k1,k5,dof*sizeof(double));

        }

        if(output) print(t,q);
        
    }

    printf("num_steps: %d\n",num_steps);
}

/** Computes a Hermite interpolation of the solution, for dense output. The
 * result is stored into the k3 array.
 * \param[in] theta the fraction of the timestep at which to evaluate the
 *                  interpolation.
 * \param[in] dt the length of the current timestep. */
void fsal_rk4d::dense_output(double theta,double dt) {
    double mth=1-theta;

    // The function assumes that the current solution is in q, the new solution
    // is in dq, the current derivative is in k2, and the new derivative is in
    // k5
    for(int i=0;i<dof;i++) k3[i]=mth*q[i]+theta*dq[i]
         -theta*mth*((1-2*theta)*(dq[i]-q[i])+dt*(theta*k5[i]-mth*k2[i]));
}

/** Performs an integration step with the fourth-order Runge-Kutta solver.
 * \param[in] dt the integration step. */
void fsal_rk4d::step(double dt) {

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+(1./3.)*dt*k1[i];
    ff(t+(1./3.)*dt,dq,k2);

    // Third RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+(-1./3.)*dt*k1[i]+dt*k2[i];
    ff(t+(2./3.)*dt,dq,k3);

    // Fourth RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*k1[i]-dt*k2[i]+dt*k3[i];
    ff(t+dt,dq,k4);
    
    // Fifth RK step - also stores derivative @ end of iterval
    for(int i=0;i<dof;i++) dq[i]=q[i]+(1./8.)*dt*(k1[i]+3.*k2[i]+3.*k3[i]+k4[i]);
    ff(t+dt,dq,k5);

    // Calculate solutions y and y^hat
    for(int i=0;i<dof;i++) {
        dq[i]=q[i]+dt*(1./8.)*(k1[i]+3.*k2[i]+3.*k3[i]+k4[i]);
        dq_hat[i]=q[i]+dt*((1./12.)*k1[i]+0.5*k2[i]+0.25*k3[i]+(1./6.)*k5[i]);
    }
    
    // Reuse k2 to store the derivative at beginning of interval
    memcpy(k2,k1,dof*sizeof(double));
    fcount+=4;
}
