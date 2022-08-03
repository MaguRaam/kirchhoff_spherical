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

    GridIn<2> gridin;
    gridin.attach_triangulation(triangulation);
    std::ifstream file("../mesh/bubble.msh");
    gridin.read_msh(file);

    Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();
    Triangulation<2>::active_cell_iterator endc = triangulation.end();

    Point<2> face_center;

    for (; cell != endc; ++cell)
        for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell; ++f)
            if (cell->face(f)->at_boundary())
                cell->face(f)->set_boundary_id(1); // Transmissive

    pcout << "===========================" << std::endl;
}
