#include "../include/Weno432.h"

Vector<double> solve_system(const FullMatrix<double>& A, const Vector<double>& b) {
    
    Vector<double> x(2);
    
    double det = A(0,0)*A(1,1) - A(0,1)*A(1,0);
    
    assert(det != 0.0);
    
    x(0) = (b(0)*A(1,1) - A(0,1)*b(1))/det;
    x(1) = (A(0,0)*b(1) - b(0)*A(1,0))/det; 
    
    return x; 
}

void Weno4_2D::output_results() {

	LA::MPI::Vector local_Ma(locally_owned_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector local_P(locally_owned_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector local_U(locally_owned_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector local_V(locally_owned_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector local_T(locally_owned_dofs, MPI_COMM_WORLD);
    LA::MPI::Vector local_Vorticity(locally_owned_dofs, MPI_COMM_WORLD);
    LA::MPI::Vector local_Schlieren(locally_owned_dofs, MPI_COMM_WORLD);

	LA::MPI::Vector Schlieren(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector Vorticity(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector Ma(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector P_out(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector U_out(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector V_out(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	LA::MPI::Vector T_out(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

	compute_primitive();
//	std::cout<<"out:\n";

//	MPI_Barrier(mpi_communicator);
//	std::cout<<"p: "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<std::endl;
    // Compute gradients using least squares approach 
       
    // Vectors and Matrices to solve the least squares gradient at cell center
    
    FullMatrix<double> A(2,2);
    Vector<double> b_RHO(2); Vector<double> b_U(2); Vector<double> b_V(2);// Vector<double> b_P(2);
    Vector<double> x_RHO(2); Vector<double> x_U(2); Vector<double> x_V(2);// Vector<double> x_P(2);
    
	double j_w, V0;
	Point<2> q_point;   
    Point<2> cell_center; 
    Point<2> neighbor_cell_center;
	unsigned int N_gp = 2, g_i, local_index; 
    
	// Weights for each neighbouring cells 

	double w, a, gamma_, p_inf_;
	double dx; 
	double dy; 
	double dRHO;
    double dU; 
    double dV;
//    double dP;
    
    // Store the summation terms 
	double w_dx_dx, w_dx_dy, w_dy_dy, w_dx_RHO, w_dy_RHO, w_dx_U, w_dy_U, w_dx_V, w_dy_V; // w_dx_P, w_dy_P
    
    Tensor<1,2> i_cap; Tensor<1,2> j_cap; 
    Tensor<1,2> r_CF; // Vector joining the centers of the two neighbouring cells 
    
    i_cap[0] = 1.0; i_cap[1] = 0.0; j_cap[0] = 0.0; j_cap[1] = 1.0;
    
	double drho_mag_old,drho_max_local = 0.0 , drho_max;
	       
    // Loop over all the cells 

	for (unsigned int c = 0; c < n_locally_cells; ++c) {

		g_i = local_to_global_index_map[c];
		local_U(g_i) = U(c);
		local_V(g_i) = V(c);
		local_P(g_i) = P(c);
		local_T(g_i) = P(c)/Rho(c);	
		gamma_ = claw.gamma(Phi(c));
		p_inf_ = claw.P_inf(Phi(c), gamma_);		
        a = std::sqrt(gamma_ * (P(c) + p_inf_) / Rho(c));
		local_Ma(g_i) = std::sqrt(U(c)*U(c) + V(c)*V(c)) / a;
//		cell = local_index_to_iterator[c];   

	        // Make the summations zero 
        
	         w_dx_dx = 0.0; w_dx_dy = 0.0; w_dy_dy = 0.0;  
    	     w_dx_RHO = 0.0; w_dy_RHO = 0.0;  
    	     w_dx_U = 0.0;  w_dy_U = 0.0;  
    	     w_dx_V = 0.0;  w_dy_V = 0.0;

            cell_center[0] = WENO_poly_consts[c](0); 
            cell_center[1] = WENO_poly_consts[c](1); 

			for (unsigned int d = 0; d < cell_all_neighbor_index[c].size(); ++d) {
//				local_neighbor_dof_indices[0] = cell_all_neighbor_index[c][d];
//				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];			
				local_index = global_to_local_index_map[cell_all_neighbor_index[c][d] ];			

				neighbor_cell_center[0] = 0.0;
				neighbor_cell_center[1] = 0.0;

				for (unsigned int i = 0; i < N_gp*N_gp; i++) {
					q_point = Cell[local_index].cell_quadrature_point(i);
					j_w = Cell[local_index].jxw(i);
					V0 = Cell[local_index].measure();
					neighbor_cell_center[0] += (1./V0)*j_w*(q_point(0));
					neighbor_cell_center[1] += (1./V0)*j_w*(q_point(1));
				}

                r_CF[0] = neighbor_cell_center(0) - cell_center(0);
				r_CF[1] = neighbor_cell_center(1) - cell_center(1);
                
                dx = r_CF*i_cap; dy = r_CF*j_cap;
                w = 1.0/neighbor_cell_center.distance(cell_center);
                
                dRHO = Rho(local_index) - Rho(c);
                dU   =   U(local_index) - U(c);
                dV   =   V(local_index) - V(c);
                
                w_dx_dx += w*dx*dx; w_dx_dy += w*dx*dy; w_dy_dy += w*dy*dy; 
                w_dx_RHO += w*dx*dRHO; w_dy_RHO += w*dy*dRHO;
                w_dx_U += w*dx*dU; w_dy_U += w*dy*dU;
                w_dx_V += w*dx*dV; w_dy_V += w*dy*dV;

			}
    
	        A(0,0) = w_dx_dx; A(0,1) = w_dx_dy; 
    	    A(1,0) = w_dx_dy; A(1,1) = w_dy_dy; 
        
    	    b_RHO(0) = w_dx_RHO; b_RHO(1) = w_dy_RHO;
    	    b_U(0)   = w_dx_U;     b_U(1) = w_dy_U;
    	    b_V(0)   = w_dx_V;     b_V(1) = w_dy_V;
        
    	    x_RHO = solve_system(A, b_RHO); x_U = solve_system(A, b_U); 
    	    x_V = solve_system(A, b_V);
    	    
			local_Schlieren(g_i) = std::sqrt( x_RHO[0]*x_RHO[0] + x_RHO[1]*x_RHO[1] );

			local_Vorticity(g_i) = x_V(0) - x_U(1);

			drho_mag_old = local_Schlieren(g_i);
			
			if( drho_mag_old > drho_max_local) 
				drho_max_local = drho_mag_old; 
    }

	drho_max = Utilities::MPI::max (drho_max_local, MPI_COMM_WORLD);

	for (unsigned int c = 0; c < n_locally_cells; ++c) {

		g_i = local_to_global_index_map[c];
		local_Schlieren(g_i) = std::exp( -1.0*local_Schlieren(g_i)/drho_max );
	}	

	local_U.compress(VectorOperation::insert);
	local_V.compress(VectorOperation::insert);
	local_P.compress(VectorOperation::insert);
	local_T.compress(VectorOperation::insert);
	local_Ma.compress(VectorOperation::insert);
	local_Vorticity.compress(VectorOperation::insert);
	local_Schlieren.compress(VectorOperation::insert);

	U_out = local_U;
	V_out = local_V;
	P_out = local_P;
	T_out = local_T;
	Ma = local_Ma;
	Vorticity = local_Vorticity;
	Schlieren = local_Schlieren;

    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);

    data_out.add_data_vector (RHO, "RHO",DataOut<2>::type_dof_data);
    data_out.add_data_vector (U_out, "U",DataOut<2>::type_dof_data);
    data_out.add_data_vector (V_out, "V",DataOut<2>::type_dof_data);
    data_out.add_data_vector (P_out, "P",DataOut<2>::type_dof_data);
    data_out.add_data_vector (Ma, "Mach",DataOut<2>::type_dof_data);
    data_out.add_data_vector (T_out, "T",DataOut<2>::type_dof_data);
    data_out.add_data_vector (PHI, "Phi",DataOut<2>::type_dof_data);
    data_out.add_data_vector (Schlieren, "Schlieren",DataOut<2>::type_dof_data);
    data_out.add_data_vector (Vorticity, "Vorticity",DataOut<2>::type_dof_data);

    data_out.build_patches ();

	unsigned int t_ = time;
	unsigned int decimal = (time - t_)*1000000;

    const std::string filename = "../plots/plot_" + Utilities::int_to_string(t_,2) + "p" + Utilities::int_to_string(decimal,6) + ".vtu";

    data_out.write_vtu_in_parallel (filename.c_str(), MPI_COMM_WORLD);
}


