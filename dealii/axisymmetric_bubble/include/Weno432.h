#ifndef WENO432_H_
#define WENO432_H_

#include "Headers.h"
#include "CLS.h"
#include "LU.h"
#include "Riemann_Solvers.h"
#include "Exceptions.h"
#include "cell_properties.h"
#include "claw.h"
#include "trueSolution.h"

// Complementary Function Definitions

const double M = 10.0;

Vector<double> initial_condition(Point<2>, double);

Vector<double> solve_system(const FullMatrix<double> &, const Vector<double> &);

double evaluate_weno_polynomial(const Vector<double> &, const Vector<double> &, const Point<2> &, double);
Point<2, double> evaluate_gradient(const Vector<double> &coeffs, const Vector<double> &consts, Point<2> P, double h);
double compute_second_order_smoothness_indicator(const Vector<double> &, const Vector<double> &, double);
double compute_third_order_smoothness_indicator(const Vector<double> &, const Vector<double> &, double);
double compute_fourth_order_smoothness_indicator(const Vector<double> &, const Vector<double> &, double);

double compute_second_order_smoothness_indicator(const Vector<double> &);
double compute_third_order_smoothness_indicator(const Vector<double> &);
double compute_fourth_order_smoothness_indicator(const Vector<double> &);

DoFHandler<2>::active_cell_iterator return_cell_pointer(const DoFHandler<2> &, unsigned int);


// Kirchhoff data:
struct KirchhoffData
{
    double theta;
    double p;
    double pr;
};

// Main Class Declaration

class Weno4_2D
{

    // functions and datastructures related to Kirchhoff surface:

    void get_quad_points();
    double get_pressure();

    // quadrature points and cells containing quadrature points:
    std::vector<Point<2, double>> q_point;
    std::vector<DoFHandler<2>::active_cell_iterator> q_cell;

    // perturbed pressure and normal derivative at quadrature points:
    std::vector<KirchhoffData> q_local, q_global;

public:
    void make_grid();
    void setup_system();
    void allocate_memory();
    void precompute_matrices();
    void precompute_matrices_veclocity();
    void compute_IS_constants();
    void compute_weno_polynomial_constants();
    void reconstruct();
    void compute_primitive();
    void initialize();
    void compute_cell_properties();
    void compute_rhs();
    void compute_time_step_based_on_cfl_number();
    void solve_ssprk33();
    void solve_ssprk54();
    void solve_ader();
    void output_results();
    void output_grid();
    void restart();
    void restart_r();
    void post_process();

private:
    MPI_Comm mpi_communicator;

    parallel::shared::Triangulation<2> triangulation;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    LA::MPI::Vector RHO;
    LA::MPI::Vector RHO_U;
    LA::MPI::Vector RHO_V;
    LA::MPI::Vector E;
    LA::MPI::Vector PHI;

    LA::MPI::Vector local_RHO;
    LA::MPI::Vector local_RHO_U;
    LA::MPI::Vector local_RHO_V;
    LA::MPI::Vector local_E;
    LA::MPI::Vector local_PHI;

    Vector<double> first_order_cell;

    Vector<double> Rho, U, V, P, Phi;

    // Coefficients for WENO polynomials
    std::vector<Vector<double>> coeffs_Rho;
    std::vector<Vector<double>> coeffs_U;
    std::vector<Vector<double>> coeffs_V;
    std::vector<Vector<double>> coeffs_P;
    std::vector<Vector<double>> coeffs_Phi;

    // WENO polynomial constants (only depend on mesh)
    std::vector<Vector<double>> WENO_poly_consts;
    std::vector<Vector<double>> IS_constants;
    std::vector<bool> is_corner_cell;
    std::vector<bool> is_relevant_cell;
    std::vector<bool> is_ghost_cell;
    std::vector<cell_properties> Cell;

    Vector<double> rhs1;
    Vector<double> rhs2;
    Vector<double> rhs3;
    Vector<double> rhs4;
    Vector<double> rhs5;

    // Fourth Order Stencil
    std::vector<Constrained_LS> CLS_R4;

    // Third Order Stencils
    std::vector<bool> is_admissible_R3; // Centered third order stencil stencil
    std::vector<Constrained_LS> CLS_R3;

    std::vector<bool> is_admissible_R31; // Third order stencil 1
    std::vector<Constrained_LS> CLS_R31;

    std::vector<bool> is_admissible_R32; // Third order stencil 2
    std::vector<Constrained_LS> CLS_R32;

    std::vector<bool> is_admissible_R33; // Third order stencil 3
    std::vector<Constrained_LS> CLS_R33;

    std::vector<bool> is_admissible_R34; // Third order stencil 4
    std::vector<Constrained_LS> CLS_R34;

    // Second order stencils

    std::vector<LUdcmp> LU_R21;
    std::vector<LUdcmp> LU_R22;
    std::vector<LUdcmp> LU_R23;
    std::vector<LUdcmp> LU_R24;

    unsigned int dofs_per_cell;
    unsigned int n_locally_cells;
    unsigned int n_vertices;
    unsigned int n_faces;
    unsigned int n_relevant_cells, n_store_cell;

    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> local_neighbor_dof_indices;

    ConditionalOStream pcout;
    TimerOutput computing_timer;

    // Slip boundaries
    std::vector<Constrained_LS> CLS_R4_slip;
    std::vector<Constrained_LS> CLS_R3_slip;

    std::map<unsigned int, unsigned int> global_to_local_index_map;
    std::map<unsigned int, unsigned int> face_index_map;
    std::vector<unsigned int> local_to_global_index_map;
    std::vector<std::vector<std::vector<unsigned int>>> cell_neighbor_index;
    std::vector<std::vector<std::vector<unsigned int>>> cell_neighbor_neighbor_index;
    std::vector<std::vector<unsigned int>> cell_all_neighbor_index;
    std::vector<std::vector<unsigned int>> cell_diagonal_neighbor_index;
    std::vector<DoFHandler<2>::active_cell_iterator> local_index_to_iterator;

    double dt;
    double finalTime;
    double cfl;
    double h_min;
    double time;

    bool RESTART, Use_ader;

    unsigned int n_refinement;

    FE_DGQ<2> fv;
    DoFHandler<2> dof_handler;
    Conservation_Law claw;
    MappingQ<2> mapping;

public:
    Weno4_2D(double, double, bool, unsigned int);
    void run();
};

#endif /* WENO432_H_ */
