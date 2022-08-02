#include "../include/Weno432.h"

void Weno4_2D::compute_rhs() {

	compute_primitive();
    reconstruct();
    
    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;

	// Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell, neighbor ;

    Point<2> face_quadrature_point_1;
    Point<2> face_quadrature_point_2;

    double V_c, h; // Volume of the cell and surface area of the face
	
	double nx1, ny1;   // Face normal vectors
	double nx2, ny2; 
	double j_w1, j_w2; 
	double gamma_L1 = 0.0, gamma_L2 = 0.0, gamma_R1 = 0.0, gamma_R2 = 0.0;
	double p_inf_L1 = 0.0, p_inf_L2 = 0.0, p_inf_R1 = 0.0, p_inf_R2 = 0.0;

	Vector<double> W(5), source(5);	
    Vector<double> WL1(5); Vector<double> WR1(5); // Solving the Riemann Problem
    Vector<double> WL2(5); Vector<double> WR2(5); // Solving the Riemann Problem
	
    Vector<double> F(5);
    
    bool boundary;

	unsigned int local_face_index;  

    std::vector< Vector<double> > Flux1(n_faces); 
    std::vector< Vector<double> > Flux2(n_faces);
    std::vector< bool > did_not_compute_flux_for_the_face(n_faces); 
    
    for (unsigned int f = 0; f < n_faces; f++) {
        Flux1[f].reinit(6);
        Flux2[f].reinit(6);
        did_not_compute_flux_for_the_face[f] = true; 
    }  
    
	unsigned int neighbor_c, g_i ;

	for (unsigned int c = 0; c < n_locally_cells; ++c) {

//		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];    

        V_c = Cell[c].measure();

        rhs1(c) = 0.0;
        rhs2(c) = 0.0;
        rhs3(c) = 0.0;
        rhs4(c) = 0.0;
        rhs5(c) = 0.0;         

		unsigned int negative = 0;
        
        for (unsigned int f = 0; f < faces_per_cell; ++f) {
            
			local_face_index = face_index_map[ cell->face_index(f) ];
			j_w1 = Cell[c].face_jxw_1(f);
			j_w2 = Cell[c].face_jxw_2(f);
            
            if(did_not_compute_flux_for_the_face[local_face_index]) {

				negative = 0;
               
                // Get some geometry info

				h = Cell[c].h(); 
    
                nx1 = Cell[c].nx1(f); ny1 = Cell[c].ny1(f);
                nx2 = Cell[c].nx2(f); ny2 = Cell[c].ny2(f);
                
                face_quadrature_point_1 = Cell[c].face_quadrature_point1(f); 
                face_quadrature_point_2 = Cell[c].face_quadrature_point2(f); 
                
                // Left face
                WL1(0) = evaluate_weno_polynomial(coeffs_Rho[c], WENO_poly_consts[c], face_quadrature_point_1, h);
                WL1(1) = evaluate_weno_polynomial(coeffs_U[c], WENO_poly_consts[c], face_quadrature_point_1, h);
                WL1(2) = evaluate_weno_polynomial(coeffs_V[c], WENO_poly_consts[c], face_quadrature_point_1, h);
                WL1(3) = evaluate_weno_polynomial(coeffs_P[c], WENO_poly_consts[c], face_quadrature_point_1, h);
                WL1(4) = evaluate_weno_polynomial(coeffs_Phi[c], WENO_poly_consts[c], face_quadrature_point_1, h);
				gamma_L1 = claw.gamma(WL1(4));
				p_inf_L1 = claw.P_inf(WL1(4), gamma_L1);

                WL2(0) = evaluate_weno_polynomial(coeffs_Rho[c], WENO_poly_consts[c], face_quadrature_point_2, h);
                WL2(1) = evaluate_weno_polynomial(coeffs_U[c], WENO_poly_consts[c], face_quadrature_point_2, h);
                WL2(2) = evaluate_weno_polynomial(coeffs_V[c], WENO_poly_consts[c], face_quadrature_point_2, h);
                WL2(3) = evaluate_weno_polynomial(coeffs_P[c], WENO_poly_consts[c], face_quadrature_point_2, h);
                WL2(4) = evaluate_weno_polynomial(coeffs_Phi[c], WENO_poly_consts[c], face_quadrature_point_2, h);

				gamma_L2 = claw.gamma(WL2(4));
				p_inf_L2 = claw.P_inf(WL2(4), gamma_L2);  
/*              
				if (WL1(0) < 0.0 || WL1(3) < 0.0 || WL2(0) < 0.0 || WL2(3) < 0.0) {
    	            std::cout<<"in compute rhs "<<std::endl<<"global index: "<<g_i<<"\t local index: "<<c
                    <<"\t rank: "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<"\tface: "<<f<<"\tquad point 1: "<<face_quadrature_point_1<<std::endl
                    <<WL1<<std::endl<<WL2<<std::endl
					<<coeffs_E[c]<<std::endl<<WENO_poly_consts[c]<<std::endl;
                }

				if (c == 0) {
    	            std::cout<<"in compute rhs "<<std::endl<<"global index: "<<c
                    <<"\tface: "<<f<<"\tquad point 1: "<<face_quadrature_point_1<<"\th: "<<h<<std::endl
                    <<WL1<<std::endl<<WL2<<std::endl
					<<coeffs_E[c]<<std::endl<<WENO_poly_consts[c]<<std::endl;
                }
*/                
                if (cell->face(f)->at_boundary()) {
                    
                    // Apply transmissive boundary conditions

	                if (cell->face(f)->boundary_id() == 0) {

						WR1 = initial_condition(face_quadrature_point_1, h_min);
						WR2 = WR1;

						gamma_R1 = claw.gamma(WR1(4));
						p_inf_R1 = claw.P_inf(WR1(4), gamma_R1);

						gamma_R2 = gamma_R1;
						p_inf_R2 = p_inf_R1;
					}
					else if (cell->face(f)->boundary_id() == 1) {

						WR1 = WL1;
						WR2 = WL2;

						gamma_R1 = claw.gamma(WR1(4));
						p_inf_R1 = claw.P_inf(WR1(4), gamma_R1);

						gamma_R2 = claw.gamma(WR2(4));
						p_inf_R2 = claw.P_inf(WR2(4), gamma_R2);

					}

					if (cell->face(f)->boundary_id() == 2) {

						WR1 = WL1;
						WR2 = WL2;

						gamma_R1 = claw.gamma(WR1(4));
						p_inf_R1 = claw.P_inf(WR1(4), gamma_R1);

						gamma_R2 = claw.gamma(WR2(4));
						p_inf_R2 = claw.P_inf(WR2(4), gamma_R2);

        	            WR1(1) = WL1(1) - 2.0*WL1(1)*nx1*nx1 - 2.0*WL1(2)*nx1*ny1; 
        	            WR1(2) = WL1(2) - 2.0*WL1(1)*nx1*ny1 - 2.0*WL1(2)*ny1*ny1;
	
        	            WR2(1) = WL2(1) - 2.0*WL2(1)*nx2*nx2 - 2.0*WL2(2)*nx2*ny2; 
        	            WR2(2) = WL2(2) - 2.0*WL2(1)*nx2*ny2 - 2.0*WL2(2)*ny2*ny2;

					}
                    
                }
                
                else {
                    
                    // Get the right state values

	    	        neighbor = cell->neighbor(f);
		            neighbor->get_dof_indices(local_neighbor_dof_indices);
					neighbor_c = global_to_local_index_map[local_neighbor_dof_indices[0] ];
					
					h = Cell[neighbor_c].h();
                    
                    WR1(0) = evaluate_weno_polynomial(coeffs_Rho[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_1, h);
                    WR1(1) = evaluate_weno_polynomial(coeffs_U[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_1, h);
                    WR1(2) = evaluate_weno_polynomial(coeffs_V[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_1, h);
                    WR1(3) = evaluate_weno_polynomial(coeffs_P[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_1, h);
                    WR1(4) = evaluate_weno_polynomial(coeffs_Phi[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_1, h);

					gamma_R1 = claw.gamma(WR1(4));
					p_inf_R1 = claw.P_inf(WR1(4), gamma_R1);
                    
                    WR2(0) = evaluate_weno_polynomial(coeffs_Rho[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_2, h);
                    WR2(1) = evaluate_weno_polynomial(coeffs_U[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_2, h);
                    WR2(2) = evaluate_weno_polynomial(coeffs_V[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_2, h);
                    WR2(3) = evaluate_weno_polynomial(coeffs_P[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_2, h);
                    WR2(4) = evaluate_weno_polynomial(coeffs_Phi[neighbor_c], WENO_poly_consts[neighbor_c],  face_quadrature_point_2, h);

					gamma_R2 = claw.gamma(WR2(4));
					p_inf_R2 = claw.P_inf(WR2(4), gamma_R2);
                
/*
		            if (WR1(0) < 0.0 || WR1(3) < 0.0 || WR2(0) < 0.0 || WR2(3) < 0.0) {
						std::cout<<"in compute rhs "<<std::endl<<"global index: "<<g_i<<"\t local index: "<<c<<std::endl<<"neighbor global index: "<<local_neighbor_dof_indices[0] <<"\t local index: "<<neighbor_c
						<<"\t rank: "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<"\tface: "<<f<<"\tcenter: "<<cell->face(f)->center()<<std::endl
						<<WR1<<std::endl<<WR2<<std::endl<<coeffs_E[neighbor_c]<<std::endl<<WENO_poly_consts[neighbor_c]<<std::endl;
					}   

					if(c == 1) {
						std::cout<<"in compute rhs "<<std::endl<<"neighbor global index: "<<neighbor_c
						<<"\tface: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\th: "<<h<<std::endl
						<<WR1<<std::endl<<WR2<<std::endl<<coeffs_E[neighbor_c]<<std::endl<<WENO_poly_consts[neighbor_c]<<std::endl;
					}                 
					if( c == 1 ) {
						std::cout<<"coeffs"<<std::endl;
						for(unsigned int i = 0; i < 10; ++i )
							std::cout<<coeffs_Rho[neighbor_c](i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 10; ++i )
							std::cout<<coeffs_Rho_U[neighbor_c](i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 10; ++i )
							std::cout<<coeffs_Rho_V[neighbor_c](i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 10; ++i )
							std::cout<<coeffs_E[neighbor_c](i)<<"\t";
						std::cout<<std::endl;
						std::cout<<"Poly"<<std::endl;
						for(unsigned int i = 0; i < 9; ++i )
							std::cout<<WENO_poly_consts[neighbor_c](i)<<"\t";	rotated_HLLC_riemann_solver(WL, WR, gamma_L, gamma_R, p_inf_L, p_inf_R, nx, ny, x, y, F_c);
						std::cout<<std::endl;							
					}		
*/                    
                    boundary = false; 
                }
/*
	            if (WL1(0) < 0.0 || (WL1(3) + p_inf_L1) < 0.0 || WR1(0) < 0.0 || (WR1(3) + p_inf_R1) < 0.0) {
					negative = 1;
				}
				unsigned int global_negative = Utilities::MPI::max (negative, MPI_COMM_WORLD);

				if(global_negative == 1) {
					pcout<<"Imaginary sound speed!\n";
					output_results();
			    }
*/
                claw.rotated_HLLC_riemann_solver(WL1, WR1, gamma_L1, gamma_R1, p_inf_L1, p_inf_R1, nx1, ny1, face_quadrature_point_1(0), face_quadrature_point_1(1), Flux1[local_face_index]);

				claw.rotated_HLLC_riemann_solver(WL2, WR2, gamma_L2, gamma_R2, p_inf_L2, p_inf_R2, nx2, ny2, face_quadrature_point_2(0), face_quadrature_point_2(1), Flux2[local_face_index]);
                    
                did_not_compute_flux_for_the_face[local_face_index] = false; 
            }
            
            else {
                
                Flux1[local_face_index] *= -1.0; 
                Flux2[local_face_index] *= -1.0;
            }
            
            
            F(0) = j_w1 * Flux1[local_face_index](0) + j_w2 * Flux2[local_face_index](0); 
   	        F(1) = j_w1 * Flux1[local_face_index](1) + j_w2 * Flux2[local_face_index](1); 
   	        F(2) = j_w1 * Flux1[local_face_index](2) + j_w2 * Flux2[local_face_index](2);
   	        F(3) = j_w1 * Flux1[local_face_index](3) + j_w2 * Flux2[local_face_index](3);
   	        F(4) = j_w1 * ( Flux1[local_face_index](4) - Phi(c)*Flux1[local_face_index](5) ) 
				 + j_w2* ( Flux2[local_face_index](4) - Phi(c)*Flux2[local_face_index](5) );

            // Add it to the rhs vectors

            rhs1(c) += (-1.0/V_c)*(F(0));
            rhs2(c) += (-1.0/V_c)*(F(1));
            rhs3(c) += (-1.0/V_c)*(F(2));
            rhs4(c) += (-1.0/V_c)*(F(3));
            rhs5(c) += (-1.0/V_c)*(F(4));
        }
        
        W(0) = Rho(c);
        W(1) = U(c);
        W(2) = V(c);
        W(3) = P(c);
        W(4) = Phi(c);
                                        
        claw.axisymmetric_source(cell->center(), W, source);

        rhs1(c) += source(0);
        rhs2(c) += source(1);
        rhs3(c) += source(2);
        rhs4(c) += source(3);
        rhs5(c) += source(4);        
        
	}
}

