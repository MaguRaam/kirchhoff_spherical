#include "../include/Weno432.h"

/*Set parameters for Kirchhoff surface (see parameters.py) */

// total no of quadarture points
int Nk = 5000;

// farfield pressure:
double p_inf = 10; // far-field pressure

// bubble radius:
double R0 = 0.25;

// radius of semicircular arc
double Rk = 25 * R0;

// observer location
double ro = 29.0 * R0;
double zo = 0.0;
Point<2, double> po{zo, ro};

// get quadrature points on Kirchhoff surface:
void Weno4_2D::get_quad_points()
{
    // store theta, coordinate of quadrature point, and cell containig quadrature points in the current process:
    double dtheta = M_PI / static_cast<double>(Nk);

    // map from theta[0, pi] -> x,y:
    double rk = Rk; // capture the globale variable:
    auto polar_to_cartesian = [rk](double theta)
    { return Point<2, double>{rk * cos(theta), rk * sin(theta)}; };

    // loop over theta and push back if the quadrature points are in current process:
    for (int i = 0; i < Nk; ++i)
    {
        double theta = (static_cast<double>(i) + 0.5) * dtheta;
        Point<2, double> p = polar_to_cartesian(theta);

        auto cell = GridTools::find_active_cell_around_point(mapping, dof_handler, p).first;

        if (cell != dof_handler.end() && cell->is_locally_owned())
        {
            q_point.push_back(p);
            q_cell.push_back(cell);
            q_local.push_back({theta, 0., 0.});
        }
    }

    // allocate memory for global data only at process 0:
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        q_global.resize(Nk);
}

// get pressure on Kirchhoff surface and return pressure at observer point:
double Weno4_2D::get_pressure()
{
    // loop over quadrature points in the current process and interpolate pressure and pressure derivative along normal:
    for (unsigned int q = 0; q < q_point.size(); ++q)
    {
        // get cell iterator:
        q_cell[q]->get_dof_indices(local_dof_indices);
        unsigned int local_index = global_to_local_index_map[local_dof_indices[0]];

        // get pressure at quadrature points:
        q_local[q].p = (evaluate_weno_polynomial(coeffs_P[local_index], WENO_poly_consts[local_index], q_point[q], Cell[local_index].h()) - p_inf);

        // gradient of pressure
        auto grad = evaluate_gradient(coeffs_P[local_index], WENO_poly_consts[local_index], q_point[q], Cell[local_index].h());

        // normal at quad point:
        Point<2, double> normal = q_point[q] / q_point[q].norm();

        // pressure derivative along normal:
        q_local[q].pr = (grad[0] * normal[0] + grad[1] * normal[1]);
    }

    // Gather local data to process 0 for writing to a file:

    //! as we are transfering struct... need a custum data type
    // 1- Here, create all the properties to call MPI_Type_create_struct
    MPI_Datatype kirchhoff_data_t;

    constexpr int ndata = 3;
    MPI_Aint displacements[ndata] = {offsetof(KirchhoffData, theta), offsetof(KirchhoffData, p), offsetof(KirchhoffData, pr)};
    int block_lengths[ndata] = {1, 1, 1};
    MPI_Datatype types[ndata] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    // 2- Create the type, and commit it
    MPI_Type_create_struct(ndata, block_lengths, displacements, types, &kirchhoff_data_t);
    MPI_Type_commit(&kirchhoff_data_t);

    // local number of sizes:
    int nproc = Utilities::MPI::n_mpi_processes(mpi_communicator);

    // size vector:
    int local_size_per_proc[nproc];
    int local_size = q_local.size();
    MPI_Allgather(&local_size, 1, MPI_INT, local_size_per_proc, 1, MPI_INT, MPI_COMM_WORLD);

    // displacement vector:
    int disp[nproc];
    int proc_begin = 0;
    for (unsigned int i = 0; i < Utilities::MPI::this_mpi_process(MPI_COMM_WORLD); ++i)
        proc_begin += local_size_per_proc[i];
    MPI_Allgather(&proc_begin, 1, MPI_INT, disp, 1, MPI_INT, MPI_COMM_WORLD);

    // mpi gather: (Not working!)
    MPI_Gatherv(q_local.data(), q_local.size(), kirchhoff_data_t, q_global.data(), local_size_per_proc, disp, kirchhoff_data_t, 0, MPI_COMM_WORLD);

    MPI_Type_free(&kirchhoff_data_t);

    // sort the vector based on theta:
    std::sort(q_global.begin(), q_global.end(), [](const KirchhoffData &a, const KirchhoffData &b)
              { return a.theta < b.theta; });

    // return pressure at observer point:
    auto cell = GridTools::find_active_cell_around_point(mapping, dof_handler, po).first;
    cell->get_dof_indices(local_dof_indices);
    unsigned int local_index = global_to_local_index_map[local_dof_indices[0]];
    return evaluate_weno_polynomial(coeffs_P[local_index], WENO_poly_consts[local_index], po, Cell[local_index].h()) - p_inf;
}
