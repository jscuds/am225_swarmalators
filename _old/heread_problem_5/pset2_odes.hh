#ifndef pset2_odes
#define pset2_odes

// Generic ODC problem class
class ode {
    public:
        /** The number of degrees of freedom. */
        int n;
        ode(int n_) : n(n_) {}
        virtual void initialize(double *y) = 0;
        virtual void equation(double t,double *y,double *dy) = 0;
};

class brusselator : public ode {
    public: 
    brusselator() : ode(2) {}
    virtual void initialize(double *y);
    virtual void equation(double t, double *y, double *dy);
};

class osc_1b : public ode {
    public: 
    osc_1b() : ode(2) {}
    virtual void initialize(double *y);
    virtual void equation(double t, double *y, double *dy);
};

class dis_7 : public ode {
    public:
        dis_7() : ode(2) {}
        virtual void initialize(double *y);
        virtual void equation(double t,double *y,double *dy);
};

class agents : public ode {
    public:
        const double J;
        const double K;
        const int N;
        agents(int N_, double J_,double K_) : N(N_), J(J_), K(K_), ode(3*N_) {}
        virtual void initialize(double *y);
        virtual void equation(double t,double *y,double *dy);
};

class agents3d : public ode {
    public:
        const double J;
        const double K;
        const int N;
        agents3d(int N_, double J_,double K_) : N(N_), J(J_), K(K_), ode(4*N_) {}
        virtual void initialize(double *y);
        virtual void equation(double t,double *y,double *dy);
};

#endif