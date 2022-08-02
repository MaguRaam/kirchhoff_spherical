#include "../include/Weno432.h"


void Weno4_2D::compute_weno_polynomial_constants() {

    pcout << "Computing WENO polynmial constants" << std::endl;

    unsigned int N_gp = 2;               // No. of quadrature points

    Point<2> q_point;

    double V0, x_0, y_0, j_w;
	
	double h, h2, h3;

	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

        V0 = Cell[c].measure();

		h = Cell[c].h();
		h2 = V0; h3 = h2*h;
		
		WENO_poly_consts[c](0) = 0.0;
		WENO_poly_consts[c](1) = 0.0; 
		WENO_poly_consts[c](2) = 0.0;
		WENO_poly_consts[c](3) = 0.0;
		WENO_poly_consts[c](4) = 0.0;
		WENO_poly_consts[c](5) = 0.0; 
		WENO_poly_consts[c](6) = 0.0;
		WENO_poly_consts[c](7) = 0.0;
		WENO_poly_consts[c](8) = 0.0;
		

		for (unsigned int i = 0; i < N_gp*N_gp; i++) {
            q_point = Cell[c].cell_quadrature_point(i);
			j_w = Cell[c].jxw(i);

			WENO_poly_consts[c](0) += (1./V0)*j_w*(q_point(0));
			WENO_poly_consts[c](1) += (1./V0)*j_w*(q_point(1));
		}
		
		x_0 = WENO_poly_consts[c](0); 
		y_0 = WENO_poly_consts[c](1);

		for (unsigned int i = 0; i < N_gp*N_gp; i++) {

            q_point = Cell[c].cell_quadrature_point(i);
			j_w = Cell[c].jxw(i);
			WENO_poly_consts[c](2) += (1./V0)*j_w*(q_point(0)-x_0)*(q_point(0)-x_0)/h2;
			WENO_poly_consts[c](3) += (1./V0)*j_w*(q_point(1)-y_0)*(q_point(1)-y_0)/h2;
			WENO_poly_consts[c](4) += (1./V0)*j_w*(q_point(0)-x_0)*(q_point(1)-y_0)/h2;
			WENO_poly_consts[c](5) += (1./V0)*j_w*(q_point(0)-x_0)*(q_point(0)-x_0)*(q_point(0)-x_0)/h3;
			WENO_poly_consts[c](6) += (1./V0)*j_w*(q_point(1)-y_0)*(q_point(1)-y_0)*(q_point(1)-y_0)/h3;
			WENO_poly_consts[c](7) += (1./V0)*j_w*(q_point(0)-x_0)*(q_point(0)-x_0)*(q_point(1)-y_0)/h3;
			WENO_poly_consts[c](8) += (1./V0)*j_w*(q_point(0)-x_0)*(q_point(1)-y_0)*(q_point(1)-y_0)/h3;
		}
    }

    pcout << "Done!" << std::endl;

}

