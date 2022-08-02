#include "../include/Weno432.h" 


// Compute constants for smoothness indicators 

void Weno4_2D::compute_IS_constants() {
    
    pcout << "Computing smoothness indicator constants" << std::endl;
    
    unsigned int N_gp = 5;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);
//    FEValues<2> fv_values (mapping, fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<2>::active_cell_iterator cell;

    Point<2> q_point;
	
	double x0, y0; 
	double x, y; 
	double h, h2, h3, h4;
    
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		cell = local_index_to_iterator[c];

		IS_constants[c] = 0.0; 

		fv_values.reinit(cell);
		h2 = 0.0; 
		for (unsigned int i=0; i<fv_values.n_quadrature_points; ++i)
            h2 += fv_values.JxW (i);	

		IS_constants[c](0) = h2;
		
		x0 = WENO_poly_consts[c](0); 
		y0 = WENO_poly_consts[c](1);
		
		h = Cell[c].h();
		h2 = h*h; h3 = h2*h; h4 = h3*h; 

		for (unsigned int i = 0; i < N_gp*N_gp; i++) {

			q_point = fv_values.quadrature_point(i); 
			
			x = q_point(0); y = q_point(1); 

			IS_constants[c](1)  += fv_values.JxW (i)*(x-x0)/h;                                 // x 
			IS_constants[c](2)  += fv_values.JxW (i)*(y-y0)/h;                                 // y
			IS_constants[c](3)  += fv_values.JxW (i)*((x-x0)*(x-x0))/h2;                       // x^2
			IS_constants[c](4)  += fv_values.JxW (i)*((y-y0)*(y-y0))/h2;                       // y^2
			IS_constants[c](5)  += fv_values.JxW (i)*((x-x0)*(y-y0))/h2;                       // xy
			IS_constants[c](6)  += fv_values.JxW (i)*((x-x0)*(x-x0)*(x-x0))/h3;                // x^3
			IS_constants[c](7)  += fv_values.JxW (i)*((y-y0)*(y-y0)*(y-y0))/h3;                // y^3 
			IS_constants[c](8)  += fv_values.JxW (i)*((x-x0)*(x-x0)*(y-y0))/h3;                // x^2y
			IS_constants[c](9)  += fv_values.JxW (i)*((x-x0)*(y-y0)*(y-y0))/h3;                // xy^2
			IS_constants[c](10) += fv_values.JxW (i)*((x-x0)*(x-x0)*(x-x0)*(x-x0))/h4;         // x^4 
			IS_constants[c](11) += fv_values.JxW (i)*((y-y0)*(y-y0)*(y-y0)*(y-y0))/h4;         // y^4
			IS_constants[c](12) += fv_values.JxW (i)*((x-x0)*(x-x0)*(y-y0)*(y-y0))/h4;         // x^2y^2
			IS_constants[c](13) += fv_values.JxW (i)*((x-x0)*(x-x0)*(x-x0)*(y-y0))/h4;         // x^3y
			IS_constants[c](14) += fv_values.JxW (i)*((x-x0)*(y-y0)*(y-y0)*(y-y0))/h4;         // xy^3 
		}
    }
    
    pcout << "Done!" << std::endl;  
}
