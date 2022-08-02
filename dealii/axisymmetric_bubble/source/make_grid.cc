#include "../include/LU.h"
#include "../include/CLS.h"
#include "../include/Headers.h"
#include "../include/Exceptions.h"
#include "../include/Riemann_Solvers.h"
#include "../include/Weno432.h"

// Make the grid

void Weno4_2D::make_grid()
{
    pcout << "Making grid" << std::endl;

    // bubble radius:
    double R = 0.25;

    // axisymmetric domain:
    double z_min = -30. * R;
    double z_max = 30. * R;
    double r_min = 0.0;
    double r_max = 30. * R;

    unsigned int n_z = 2000;
    unsigned int n_r = 1000;
    std::vector<unsigned int> repetions{n_z, n_r};

    bool colorize = false; // Set boundary ids for the four boundaries

    Point<2> P1(z_min, r_min);
    Point<2> P2(z_max, r_max);
    GridGenerator::subdivided_hyper_rectangle(triangulation, repetions, P1, P2, colorize);

    Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();
    Triangulation<2>::active_cell_iterator endc = triangulation.end();

    Point<2> face_center;

    for (; cell != endc; ++cell)
        for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
            if (cell->face(f)->at_boundary())
                cell->face(f)->set_boundary_id(1); // Transmissive

    pcout << "===========================" << std::endl;
}
