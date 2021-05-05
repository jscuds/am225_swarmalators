#ifndef FSAL_RK4D_HH
#define FSAL_RK4D_HH

#include<string>

/** Class for solving an ODE IVP using a 4th order FSAL method (built on rk4 solver). */
class fsal_rk4d {
    public:
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The current time. */
        double t;
        /** The solution vector. */
        double *q;
        /** Error tolerance for step size*/
        double Atol;

        std::string filename = "test.txt";

        fsal_rk4d(int dof_);
        virtual ~fsal_rk4d();
        void print(double t_,double *in);
        void print_state(double t_,double *in);
        void dense_output(double theta,double dt);
        void solve_fixed(double t_end,double Atol,bool output,int d_steps);
        void step(double dt);
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    private:
        double *dq;
        double *dq_hat;
        double *k1;
        double *k2;
        double *k3;
        double *k4;
        double *k5;
};

#endif
