/*
 * hype1d.cc
 *      Author: sunder
 */

//----------------------------------------------------------------------------
// Commonly used C/C++ header files
//----------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <ctime>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <array>
#include <map>
#include <utility>
//#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

using namespace boost; 

//----------------------------------------------------------------------------
// Various constants used throught the code
//----------------------------------------------------------------------------

const int nVar = 4;             // Number of components in the PDE system //
const int nLin = 2;             // Number of linear degenerate fields in the PDE system //

const double gamma_1       = 4.4;     // Specific heat ratio of the liquid phase //
const double gamma_2       = 1.4;     // Specific heat ratio of the gas phase //
const double P_1           = 600.0;   // Stiffness constant of the liquid phase //
const double P_2           = 0.0;     // Stiffness constant of the gas phase //
const double prs_floor     = 1.0e-12; // Pressure floor value //
const double rho_floor     = 1.0e-14; // Density floor value //
const double small_num     = 1.0e-12; // Effective small number in the code //
const int  dofs_per_cell   = 2;       // Number of degrees of freedom polynomial expansion in a cell //

const double R0            = 1.0;   // Bubble radius //

typedef multi_array<double, 2> Matrix;
typedef multi_array<double, 1> Vector;

//----------------------------------------------------------------------------
// Various types of boundary conditions
//----------------------------------------------------------------------------

enum bndry_type{inflow, periodic, reflective, transmissive};

//----------------------------------------------------------------------------
// Structure defining various critical parameters controlling the simulation
//----------------------------------------------------------------------------

struct  AppCtx {
    double x_min;                      /* x-coordinate of the domain begining */
    double x_max;                      /* x-coordinate of the domain ending */
    double N_cells;                    /* No. of cells in the domain */
    double CFL;                        /* CFL condition, should be less than 1.0 */
    double InitialTime = 0.0;          /* Initial time of the simulation */
    double FinalTime;                  /* Final time of the simulation */
    int write_interval;                /* Number of time steps aftter which to write output file */
    bool reconstruct_primitive = true; /* Reconstruct primitive variables */

    enum bndry_type left_boundary;     /* Boundary condition on the left face */
    enum bndry_type right_boundary;    /* Boundary condition on the right face */
};

/* -------------------------------------------------------- Linear-Algebra Functions -------------------------------------------------------- */

//----------------------------------------------------------------------------
// Find the minimum value in a vector V of size N 
//----------------------------------------------------------------------------

double MinVal(const Vector& V, const int& N) {

    double min = V[0]; 
    
    for (int i = 1; i < N; ++i) {
        if (V[i] < min)
            min = V[i]; 
    }
    
    return min; 
}

//----------------------------------------------------------------------------
// Find the maximum value in a vector V of size N 
//---------------------------------------------------------------------------

double MaxVal(const Vector& V, const int& N) {

    double max = V[0]; 
    
    for (int i = 1; i < N; ++i) {
        if (V[i] > max)
            max = V[i]; 
    }
    
    return max; 
}

//----------------------------------------------------------------------------
// Given two vectors V1 and V2 of size N, find the L2 residue between them  
//---------------------------------------------------------------------------

double residue(const Vector& V1, const Vector& V2, const int& N) {
    
    double sum = 0.0;
    double diff; 
    
    for (int i = 0; i < N; ++i) {
        diff = V1[i] - V2[i];
        sum += diff*diff; 
    }
    
    return std::sqrt(sum); 
}

//----------------------------------------------------------------------------
// Matrix-Vector Multiplication
// Does the operation: y := A*x
// A -> (m x n) matrix
// x -> (n) vector
// y -> (m) vector
//----------------------------------------------------------------------------

void MatVecMult(const Matrix& A, const Vector& x, Vector& y, int m, int n) {
    for (int i = 0; i < m; i++ ) {
        y[i]= 0.0;
        for (int j = 0; j < n; j++ )
            y[i]+= A[i][j]*x[j];
    }
}

//----------------------------------------------------------------------------
// Matrix-Matrix Multiplication
// Does the operation: C := A*B
// A -> (n x m) matrix
// B -> (m x p) matrix
// C -> (n x p) matrix
//----------------------------------------------------------------------------

void MatMatMult(const Matrix& A, const Matrix& B, Matrix& C, int n, int m, int p) {
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            C[i][j] = 0.0; 
            for (int k = 0; k < m; ++k) {
                C[i][j] += A[i][k]*B[k][j]; 
            }
        }
    }
}

/* -------------------------------------------------------- PDE Functions -------------------------------------------------------- */

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable
//----------------------------------------------------------------------------

void Cons2Prim(const Vector& Q, Vector& V) {

    double g = 1.0 + (gamma_1-1.0)*(gamma_2-1.0)/((1.0-Q[3])*(gamma_1 - 1.0) + Q[3]*(gamma_2-1.0));
    double P_inf = ((g -1.0)/g)*( gamma_1*P_1*Q[3]/(gamma_1 - 1.0) + gamma_2*P_2*(1 - Q[3])/(gamma_2 - 1.0) );
    
    V[0] = Q[0];
    V[1] = Q[1]/Q[0];
    V[2] = (g -1.0)*( Q[2] - 0.5*V[0]*V[1]*V[1])   - g*P_inf;
    V[3] = Q[3];  

}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
//----------------------------------------------------------------------------

void Prim2Cons(const Vector& V, Vector& Q) {

    double g = 1.0 + (gamma_1-1.0)*(gamma_2-1.0)/((1.0-V[3])*(gamma_1 - 1.0) + V[3]*(gamma_2-1.0));
    double P_inf = ((g -1.0)/g)*( gamma_1*P_1*V[3]/(gamma_1 - 1.0) + gamma_2*P_2*(1 - V[3])/(gamma_2 - 1.0) );

    Q[0] = V[0]; 
    Q[1] = V[0]*V[1];
    double e = (V[2] + g*P_inf)/(g - 1.0);  
    double k = 0.5*V[0]*V[1]*V[1];
    Q[2] = k + e; 
    Q[3] = V[3]; 

}

//----------------------------------------------------------------------------
// Find the conservative flux components F in the given normal direction
//----------------------------------------------------------------------------

double PDEFlux(const Vector& Q, const double& x, Vector& F) {

    double g = 1.0 + (gamma_1-1.0)*(gamma_2-1.0)/((1.0-Q[3])*(gamma_1 - 1.0) + Q[3]*(gamma_2-1.0));
    double P_inf = ((g -1.0)/g)*( gamma_1*P_1*Q[3]/(gamma_1 - 1.0) + gamma_2*P_2*(1 - Q[3])/(gamma_2 - 1.0) );
    
    double rho = Q[0];
    double u = Q[1]/Q[0];
    double p = (g -1.0)*( Q[2] - 0.5*rho*u*u)   - g*P_inf;

    // Check if the input state is physically admissible

    if (rho < rho_floor) {
        std::cerr << "Negative density = " << rho << std::endl;
        std::cerr << "At x = " << x << std::endl;
        std::exit(EXIT_FAILURE);
    }


    if ((p + P_inf)  < prs_floor) {
        std::cerr << "Negative pressure = " << p + P_inf << std::endl;
        std::cerr << "At x = " << x << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Now find the fluxes

    F[0] = rho*u;
    F[1] = rho*u*u + p;
    F[2] = u*(Q[2] + p);
    F[3] = 0.0;

    // Return maximum eigen value
    
    return (std::abs(u) + std::sqrt(g*(p + P_inf)/rho));
}

//----------------------------------------------------------------------------
// Find the non-conservative flux components nF
//----------------------------------------------------------------------------

void PDENonConsFlux(const Vector& Q, const Vector& grad_Q_x, Vector& nF) {
    
    double u = Q[1]/Q[0]; 
    double phi_x = grad_Q_x[3];

    // Now find the fluxes

    nF[0] = 0.0;
    nF[1] = 0.0;
    nF[2] = 0.0;
    nF[3] = u*phi_x;
}

//----------------------------------------------------------------------------
// For a given state Q, find all the eigenvalues L
//----------------------------------------------------------------------------

void PDEEigenvalues(const Vector& Q, Vector& L) {
    
    double g = 1.0 + (gamma_1-1.0)*(gamma_2-1.0)/((1.0-Q[3])*(gamma_1 - 1.0) + Q[3]*(gamma_2-1.0));
    double P_inf = ((g -1.0)/g)*( gamma_1*P_1*Q[3]/(gamma_1 - 1.0) + gamma_2*P_2*(1 - Q[3])/(gamma_2 - 1.0) );
    
    double rho = Q[0];
    double u = Q[1]/Q[0];
    double p = (g -1.0)*( Q[2] - 0.5*rho*u*u)   - g*P_inf;
    
    double a = std::sqrt(g*(p + P_inf)/rho);
    
    L[0] = u;       
    L[1] = u;               
    L[2] = u - a;        
    L[3] = u + a; 
}


//----------------------------------------------------------------------------
// PDE Source Term
//----------------------------------------------------------------------------

void PDESource(double x, const Vector& Q, Vector& S) {

    double g = 1.0 + (gamma_1-1.0)*(gamma_2-1.0)/((1.0-Q[3])*(gamma_1 - 1.0) + Q[3]*(gamma_2-1.0));
    double P_inf = ((g -1.0)/g)*( gamma_1*P_1*Q[3]/(gamma_1 - 1.0) + gamma_2*P_2*(1 - Q[3])/(gamma_2 - 1.0) );

    double rho = Q[0];
    double u = Q[1]/Q[0];
    double E = Q[2];
    double p = (g -1.0)*( E - 0.5*rho*u*u)   - g*P_inf;


    S[0] = -2.0*rho*u/x;
    S[1] = S[0]*u;
    S[2] = -2.0*(E+p)*u/x;
    S[3] = 0.0;

}

//----------------------------------------------------------------------------
// Form the non-conservative matrix B
//----------------------------------------------------------------------------

void PDEMatrixB(const Vector& Q, Matrix& B) {

    for (int i = 0; i < nVar; ++i)
        for (int j = 0; j < nVar; ++j)
            B[i][j] = 0.0;

    B[3][3] =  Q[1]/Q[0]; // u
}

//----------------------------------------------------------------------------
// The Roe matrix B̃(Ua,Ub) between two generic states Qa and Qb is numerically
// via Gaussian quadrature, by the function RoeMatrix
//----------------------------------------------------------------------------

void RoeMatrix(const Vector& Qa, const Vector& Qb, Matrix& BRoe) {

    // 3 - point quadrature points in [0,1] for path integrals

    const int N_gps = 3;
    const double s_gp[] = {0.1127016653792583, 0.5, 0.8872983346207417};
    const double w_gp[] = {0.2777777777777778, 0.4444444444444444, 0.2777777777777778};
    
    Matrix B(extents[nVar][nVar]);
    Vector Q(extents[nVar]);

    // First make the BRoe matrix zero

    for (int i = 0; i < nVar; ++i)
        for (int j = 0; j < nVar; ++j)
            BRoe[i][j] = 0.0;

    for (int q = 0; q < N_gps; ++q) {

        for (int c = 0; c < nVar; ++c)
            Q[c] = Qa[c] + s_gp[q]*(Qb[c] - Qa[c]);

        PDEMatrixB(Q, B);

        for (int i = 0; i < nVar; ++i) {
            for (int j = 0; j < nVar; ++j) {
                BRoe[i][j] += w_gp[q]*B[i][j];
            }
        }
    }
}

//----------------------------------------------------------------------------
// LLF (Rusanov) Non-Conservative Flux 
//----------------------------------------------------------------------------

double LLFNC(const Vector& QL, const Vector& QR, const double& x, Vector& F, Vector& D) {

    Vector FL(extents[nVar]), FR(extents[nVar]), Q_jump(extents[nVar]);
    Matrix B(extents[nVar][nVar]);

    double s_max_l = PDEFlux(QL, x, FL);
    double s_max_r = PDEFlux(QR, x, FR);

    double s_max = std::max(s_max_l, s_max_r);

    for (int c = 0; c < nVar; ++c) {
        Q_jump[c] = QR[c] - QL[c];
        F[c] = 0.5*(FR[c] + FL[c] - s_max*(Q_jump[c]));
    }

    RoeMatrix(QL, QR, B);
    
    // Multiply B with the jump term

    MatVecMult(B, Q_jump, D, nVar, nVar);

    return s_max;
}

//----------------------------------------------------------------------------
// Minmod limiter function 
//----------------------------------------------------------------------------

double sgn(double a) {
    if (a > 0.0) return 1.0;
    else if ( a == 0.0) return 0.0;
    else return -1.0;
}

double minmod_limiter(double a, double b) {
    return 0.5*(sgn(a)+ sgn(b))*std::min(std::abs(a), std::abs(b));
}

//----------------------------------------------------------------------------
// Initial Condition Function
//----------------------------------------------------------------------------


double TanhRadial(double R, double R0, double eps, double in, double out)
{
    return 0.5 * ((in + out) + (out - in) * tanh((R - R0) / (eps)));
}

Vector initial_condition(double x) {

    Vector V0(extents[nVar]);

    double smear = 1.0;
    double h = .01;

    double p_medium = 1.0;  // pressure in the water medium
    double p_b = 0.1;       // Pressure in the air bubble

    double rho_medium = 1.0; // density of the water medium
    double rho_b = 0.001;    // density of the air bubble

    V0[0] = TanhRadial(x, R0, smear*h, rho_b, rho_medium);
    V0[1] = 0.0; 
    V0[2] = TanhRadial(x, R0, smear*h, p_b, p_medium);
    V0[3] = TanhRadial(x, R0, smear*h, 0.0, 1.0);

    return V0; 
}

//----------------------------------------------------------------------------
// Path Conservative Finite Volume Method Class
//----------------------------------------------------------------------------

class HyPE_1D {

    // Typedefs

    typedef boost::multi_array<Vector, 1> array_type;
    typedef array_type::index index;

    boost::multi_array<std::vector<double>, 2> U;  // Conserved variables at cells 
    boost::multi_array<Vector, 1> Dp;              // Non-Conservative fluxes at cell faces
    boost::multi_array<Vector, 1> Dm;              // Non-Conservative fluxes at cell faces
    boost::multi_array<Vector, 1> RHS;             // RHS term for each face 
    boost::multi_array<Vector, 1> U_L;             // Value of conserved variable at cell left face 
    boost::multi_array<Vector, 1> U_R;             // Value of conserved variable at cell right face 
    boost::multi_array<Vector, 1> F_L;             // Value of conserved variable at cell left face 
    boost::multi_array<Vector, 1> F_R;             // Value of conserved variable at cell right face 
    boost::multi_array<double, 1> x;               // Cell centers

    AppCtx Params;
    int N_ph;
    double dx;
    double dt;
    double time;
    int time_step;
    int rk_stage;

    void initialize();
    void apply_boundary_conditions();
    void limit_solution();
    void compute_rhs(double);
    void solve();
    void plot(int, unsigned int = 4) const;
    void compute_bubble_radius() const;

public:
    HyPE_1D(AppCtx);
    void run();
};

//----------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------

HyPE_1D::HyPE_1D(AppCtx params) :
    U(boost::extents[params.N_cells + 6][nVar]),
    Dp(boost::extents[params.N_cells + 1]),
    Dm(boost::extents[params.N_cells + 1]),
    RHS(boost::extents[params.N_cells]),
    U_L(boost::extents[params.N_cells + 2]),
    U_R(boost::extents[params.N_cells + 2]),
    F_L(boost::extents[params.N_cells + 2]),
    F_R(boost::extents[params.N_cells + 2]),
    x(boost::extents[params.N_cells + 6]),
    Params(params),
    N_ph(3),
    dt(0.0),
    time(0.0),
    time_step(0),
    rk_stage(0)
    {

    // Initialize the grid

    dx = (Params.x_max - Params.x_min)/static_cast<double>(Params.N_cells);
    boost::array<array_type::index, 1> bases_1d = {{-N_ph}};
    x.reindex(bases_1d);

    for (int i = -N_ph; i < Params.N_cells + N_ph; i++)
        x[i] = Params.x_min + (static_cast<double>(i)+0.5)*dx;

    boost::array<array_type::index, 2> bases_2d = {{-N_ph, 0}};

    U.reindex(bases_2d);


    for (int i = -3; i < Params.N_cells + 3; ++i)
        for (int c = 0; c < nVar; ++c)
                U[i][c].resize(dofs_per_cell, 0.0);

    boost::array<array_type::index, 1> bases_1d_face = {{-1}};

    U_L.reindex(bases_1d_face);
    U_R.reindex(bases_1d_face);
    F_L.reindex(bases_1d_face);
    F_R.reindex(bases_1d_face);
    
    for (int i = 0; i < Params.N_cells; ++i) {
        
        RHS[i].resize(extents[nVar]);
    }
    
    for (int i = -1; i < Params.N_cells+1; ++i) {
        
        U_L[i].resize(extents[nVar]);
        U_R[i].resize(extents[nVar]);
        F_L[i].resize(extents[nVar]);
        F_R[i].resize(extents[nVar]);
        
    }
    
    for (int i = 0; i < Params.N_cells + 1; ++i) {
            
        Dp[i].resize(extents[nVar]);
        Dm[i].resize(extents[nVar]);
        
    }
}

//----------------------------------------------------------------------------
// Initialize the solution using midpoint rule 
//----------------------------------------------------------------------------

void HyPE_1D::initialize() {
    
    std::cout << "Initializing the solution" << std::endl;

    Vector Q(extents[nVar]), V(extents[nVar]);

    // Loop through all the cells

    for (int i = 0; i < Params.N_cells; i++) {

        V = initial_condition(x[i]);

        Prim2Cons(V, Q);
        
        for (int c = 0; c < nVar; ++c)
            U[i][c][0] = Q[c];
    }

    std::cout << "Done!" << std::endl;
}

//----------------------------------------------------------------------------
// Apply boundary conditions
//----------------------------------------------------------------------------

void HyPE_1D::apply_boundary_conditions() {

    int oned_begin, oned_end, ilhs, irhs;

    // ---------------------- Left boundary ----------------------

    oned_begin = 0; oned_end = Params.N_cells-1;

    for (int i = 0; i < N_ph; ++i) {

        // Outflow/Transmissive boundary

        if (Params.left_boundary == transmissive) {

            ilhs = oned_begin - i - 1;
            irhs = oned_begin;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];
            
        }

        if (Params.left_boundary == reflective) {

            ilhs = oned_begin - i - 1;
            irhs = oned_begin + i;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];

            U[ilhs][1][0] = -U[ilhs][1][0];
            U[ilhs][4][0] = -U[ilhs][4][0];
        }

        if (Params.left_boundary == periodic) {
            ilhs = oned_begin - i - 1;
            irhs = oned_end - i;

            for (int c = 0; c < nVar; ++c)
                U[ilhs][c][0] = U[irhs][c][0];
            
        }
    }

    // ---------------------- Right boundary ----------------------

    for (int i = 0; i < N_ph; ++i) {

        // Outflow/Transmissive boundary

        if (Params.right_boundary == transmissive) {

            ilhs = oned_end + i + 1;
            irhs = oned_end;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];
        }

        if (Params.right_boundary == reflective) {

            ilhs = oned_end + i + 1;
            irhs = oned_end - i;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];

            U[ilhs][1][0] = -U[ilhs][1][0];
        }

        if (Params.right_boundary == periodic) {

            ilhs = oned_end + i + 1;
            irhs = oned_begin + i;

            for (int c = 0; c < nVar; ++c) 
                U[ilhs][c][0] = U[irhs][c][0];
        }
    }
}

//----------------------------------------------------------------------------
// Find the RHS in each cell
//----------------------------------------------------------------------------

void HyPE_1D::compute_rhs(double t) {

    double s, left_slope, right_slope, s_max = 0.0;
    const double r1_dx = 1./dx;
    Vector Q(extents[nVar]);
    Vector grad_Q(extents[nVar]);
    Vector nF(extents[nVar]);
    Vector Source(extents[nVar]);
    Vector Q_L(extents[nVar]); Vector Q_R(extents[nVar]);
    
    apply_boundary_conditions();

    for (int i = -1; i < Params.N_cells+1; ++i) {

        for (int c = 0; c < nVar; ++c) {

            left_slope  = U[i][c][0] - U[i-1][c][0];
            right_slope = U[i+1][c][0] - U[i][c][0];
            
            U[i][c][1] = minmod_limiter(left_slope, right_slope);
            
            U_L[i][c] = U[i][c][0] - 0.5*U[i][c][1];
            U_R[i][c] = U[i][c][0] + 0.5*U[i][c][1];
        }
        
        PDEFlux(U_L[i], x[i]-0.5*dx, F_L[i]); 
        PDEFlux(U_R[i], x[i]+0.5*dx, F_R[i]); 
    }

    // Find upwind flux

    for (int i = 0; i < Params.N_cells + 1; ++i) {

        for (int c = 0; c < nVar; ++c) {
            Q_L[c] = U_R[i-1][c];
            Q_R[c] = U_L[i][c];
        }
    
        s = LLFNC(Q_L, Q_R, x[i]-0.5*dx, Dm[i], Dp[i]); 

        if (s > s_max)
            s_max = s;
    }

    // Find RHS

    for (int i = 0; i < Params.N_cells; ++i) 
        for (int c = 0; c < nVar; ++c) 
            RHS[i][c] = -r1_dx*(Dm[i+1][c] - Dm[i][c] + 0.5*(Dp[i+1][c] + Dp[i][c]));
    
    // Add the smooth part of the non-conservative term

    for (int i = 0; i < Params.N_cells; ++i) {

        for (int c = 0; c < nVar; ++c) {
            Q[c] = U[i][c][0];
            grad_Q[c] = r1_dx*U[i][c][1]; 
        }

        PDENonConsFlux(Q, grad_Q, nF);
        PDESource(x[i],Q,Source);

        for (int c = 0; c < nVar; ++c)
            RHS[i][c] += (Source[c] - nF[c]);
    }

    // Find the time step size (only if rk_stage = 1)

    if (rk_stage == 1) {
        dt = Params.CFL*dx/s_max;

        // If time step exceeds the final time, reduce it accordingly

        if((time + dt)>Params.FinalTime)
            dt = Params.FinalTime - time;
    }
}

//----------------------------------------------------------------------------
// Update Solution using SSPRK22 method
//----------------------------------------------------------------------------

void HyPE_1D::solve() {

    std::cout << "Solving using SSPRK (2,2) method" << std::endl;

    boost::multi_array<double, 2> U_old(extents[Params.N_cells][nVar]);

    while (time < Params.FinalTime) {

        // write time and bubble radius:
        compute_bubble_radius();

        printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.FinalTime);

        //  Stage 1

        rk_stage = 1;
        compute_rhs(time);

        for (int i = 0; i < Params.N_cells; ++i){

            for (int c = 0; c < nVar; ++c) {

                U_old[i][c] =   U[i][c][0];

                U[i][c][0] +=   dt*RHS[i][c];

            }
        }
        
        //  Stage 2

        rk_stage = 2;
        compute_rhs(time + dt);


        for (int i = 0; i < Params.N_cells; ++i) {

            for (int c = 0; c < nVar; ++c)
                U[i][c][0] =   0.5*(U_old[i][c] +   U[i][c][0] + dt*RHS[i][c]);

        }

        time += dt;
        time_step++;
    }

    printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.FinalTime);
}

//----------------------------------------------------------------------------
// Plot solution as csv file
//----------------------------------------------------------------------------

std::string int_to_string (int value, const unsigned int digits) {
    std::string lc_string = std::to_string(value);

    if (lc_string.size() < digits) {
        // We have to add the padding zeroes in front of the number
        const unsigned int padding_position = (lc_string[0] == '-')
                                                ?
                                                1
                                                :
                                                0;

        const std::string padding(digits - lc_string.size(), '0');
        lc_string.insert(padding_position, padding);
    }

    return lc_string;
}


void HyPE_1D::plot(int i, unsigned int digits) const {

    Vector V(extents[nVar]), Q(extents[nVar]);

    std::ofstream out_data;
    const std::string filename = "plot/sol-" + int_to_string (i, digits) + ".csv";
    out_data.open (filename);
    out_data.flags( std::ios::dec | std::ios::scientific );
    out_data.precision(6);

    //out_data << "x,RHO,V,P,PHI" << std::endl;

    for (int i = 0; i < Params.N_cells; ++i) {

        for (int c = 0; c < nVar; ++c) {
            Q[c] = U[i][c][0];
        }

        Cons2Prim(Q, V);

        out_data << x[i] << ",";

        for (int c = 0; c < nVar; ++c) {

            if (c == nVar -1)
                out_data << V[c];

            else
                out_data << V[c] << ",";
        }

        out_data << std::endl;

    }

    out_data.close();
}


//----------------------------------------------------------------------------
// Compute bubble radius
//----------------------------------------------------------------------------
void HyPE_1D::compute_bubble_radius() const {

    std::ofstream file("bubble_radius.csv", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    double bubble_radius = 0.0;

    for (int i = 0; i < Params.N_cells; ++i)
        bubble_radius += (1.0 - U[i][3][0])*dx;

    file << time << "," << bubble_radius << "\n";
    
}


//----------------------------------------------------------------------------
// Put everything together and run the problem
//----------------------------------------------------------------------------

void HyPE_1D::run() {
    
    auto start = std::chrono::system_clock::now();

    //-------------------------------------
    
    initialize();
    solve();
    plot(time_step);
    
    //-------------------------------------
    
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "No. of time steps = " << time_step << std::endl; 
    std::cout << "Time taken = " << elapsed_seconds.count() << std::endl;
}

//----------------------------------------------------------------------------
// Main function of the code 
//----------------------------------------------------------------------------

int main() {

    AppCtx Params;

    Params.x_min = 0.0;
    Params.x_max = 10*R0;
    Params.CFL   = 0.8;
    Params.InitialTime = 0.0;
    Params.FinalTime = 1.68;
    Params.N_cells = 1000;
    Params.write_interval = 40;
    Params.left_boundary  = transmissive;
    Params.right_boundary = transmissive;

    HyPE_1D Effective_Gamma(Params);

    Effective_Gamma.run();
    
    return 0;
}
