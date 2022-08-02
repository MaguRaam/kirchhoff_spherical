#include "../include/Weno432.h"

// Precompute reconstruction matrices 

void Weno4_2D::precompute_matrices() {

    TimerOutput::Scope t(computing_timer, "Computing the stencil matrices");

//    auto start = std::chrono::system_clock::now();

    pcout << "Computing the stencil matrices" << std::endl;

    unsigned int N_gp = 2;  // No. of quadrature points

    Point<2> q_point;
    Point<2> C; 
    Point<2> P; 
    
    double V_neighbor; 

    DoFHandler<2>::active_cell_iterator cell, neighbor;


    FullMatrix<double> A_R4; // Fourth Order Stencil (Least Squares Part)
    FullMatrix<double> C_R4; // Fourth Order Stencil (Constraint Part)
    
    FullMatrix<double> A_R3; // Least Squares Matrix for r=3 stencil 
    FullMatrix<double> C_R3; // Constraint Matrix for r=3 stencil 
    
    // One-sided stencils 
    
	FullMatrix<double> C_R31; 
    FullMatrix<double> A_R31;
    
    FullMatrix<double> C_R32; 
    FullMatrix<double> A_R32;
    
    FullMatrix<double> C_R33; 
    FullMatrix<double> A_R33;

	FullMatrix<double> C_R34; 
    FullMatrix<double> A_R34;
    
	FullMatrix<double> A_R21(2,2); // Matrix for r=2 stencil 1
    FullMatrix<double> A_R22(2,2); // Matrix for r=2 stencil 2
    FullMatrix<double> A_R23(2,2); // Matrix for r=2 stencil 3
    FullMatrix<double> A_R24(2,2); // Matrix for r=2 stencil 4
    
    unsigned int ROWS, index, WW_local_index, NN_local_index, EE_local_index, SS_local_index, local_index;
	
	double x0, y0, h, h2, h3, j_w; 
    
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		cell = local_index_to_iterator[c];

		bool WW = false, NN = false, SS = false, EE = false;
		
		x0 = WENO_poly_consts[c](0); 
        y0 = WENO_poly_consts[c](1);
		h = Cell[c].h();
		h2 = h*h; h3 = h2*h;

	
		if (cell_neighbor_neighbor_index[c][0].size() > 0) {
			WW = true;
			WW_local_index = global_to_local_index_map[cell_neighbor_neighbor_index[c][0][0] ];
		}

		if (cell_neighbor_neighbor_index[c][1].size() > 0) {
			EE = true;
			EE_local_index = global_to_local_index_map[cell_neighbor_neighbor_index[c][1][0] ];
		}

		if (cell_neighbor_neighbor_index[c][2].size() > 0) {
			SS = true;
			SS_local_index = global_to_local_index_map[cell_neighbor_neighbor_index[c][2][0] ];
		}

		if (cell_neighbor_neighbor_index[c][3].size() > 0) {
			NN = true;
			NN_local_index = global_to_local_index_map[cell_neighbor_neighbor_index[c][3][0] ];
		}
       
        if ( !(cell->at_boundary()) ) {
			
			// =====================================================================
            // r = 4 stencil 
            // =====================================================================
			
			C_R4.reinit(4, 9); C_R4 = 0.0;
            
            x0 = WENO_poly_consts[c](0); 
            y0 = WENO_poly_consts[c](1);
            
            // Fill C_R4 (Constraint r = 4 matrix)

            // Row 1 (W neighbor)
            
			neighbor = cell->neighbor(0);
		    neighbor->get_dof_indices(local_neighbor_dof_indices);
			local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	V_neighbor = Cell[local_index].measure();
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                C_R4(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                C_R4(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                C_R4(0,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                C_R4(0,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                C_R4(0,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                C_R4(0,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                C_R4(0,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                C_R4(0,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                C_R4(0,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
            }

            // Row 2 (N neighbor)
			neighbor = cell->neighbor(3);
            neighbor->get_dof_indices(local_neighbor_dof_indices);
			local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
			V_neighbor = Cell[local_index].measure();
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                C_R4(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                C_R4(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                C_R4(1,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                C_R4(1,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                C_R4(1,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                C_R4(1,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                C_R4(1,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                C_R4(1,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                C_R4(1,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
            }

            // Row 3 (E neighbor) 
            
			neighbor = cell->neighbor(1);
            neighbor->get_dof_indices(local_neighbor_dof_indices);
			local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
			V_neighbor = Cell[local_index].measure();
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                C_R4(2,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                C_R4(2,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                C_R4(2,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                C_R4(2,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                C_R4(2,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                C_R4(2,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                C_R4(2,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                C_R4(2,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                C_R4(2,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
            }

            // Row 4 (S neighbor) 
            
			neighbor = cell->neighbor(2);
            neighbor->get_dof_indices(local_neighbor_dof_indices);
			local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
			V_neighbor = Cell[local_index].measure();
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                C_R4(3,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                C_R4(3,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                C_R4(3,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                C_R4(3,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                C_R4(3,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                C_R4(3,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                C_R4(3,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                C_R4(3,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                C_R4(3,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
            }
            
            // Least Squares Part  
            
			// First get the number of rows in the least squares matrix part 
			
			// vertex neighbors 

			ROWS = cell_diagonal_neighbor_index[c].size(); 
			
			// neighbor of neighbors 

			if (WW) {
				ROWS++; 
			}

			if (EE) {
				ROWS++; 
			}

			if (SS) {
				ROWS++; 
			}

			if (NN) {
				ROWS++; 
			}
				
			index = 0; 
            
            A_R4.reinit(ROWS, 9); A_R4 = 0.0;
		
			// First fill the vertex neighbors 
			
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			
	
				local_index = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];
    	        V_neighbor = Cell[local_index].measure();
				
				for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                    A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                    A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                    A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
                }
                
                index++; 
			}
            
			if (WW) {
				
				local_index = WW_local_index;
	            V_neighbor = Cell[local_index].measure();
				
				for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                    A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                    A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                    A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
                }
                
                index++; 
			}
			
			if (NN) {
				
				local_index = NN_local_index;
	            V_neighbor = Cell[local_index].measure();
				
				for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                    A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                    A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                    A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
                }
                
                index++; 
			}
			
			if (EE) {
				
				local_index = EE_local_index;
	            V_neighbor = Cell[local_index].measure();
				
				for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                    A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                    A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                    A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
                }
                
                index++; 
			}
			
			if (SS) {
				
				local_index = SS_local_index;
	            V_neighbor = Cell[local_index].measure();
				
				for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
                    A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
                    A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
                    A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
                }
                
                index++; 
			}
            
            CLS_R4[c].initialize(A_R4, C_R4); 
            
            // =====================================================================
            // r = 3 stencil (Centered Stencil)
            // =====================================================================
            
			C_R3.reinit(4, 5); C_R3 = 0.0;
			
			// Fill C_R3 (Constraint r = 3 matrix)

            // Row 1 (W neighbour)

			neighbor = cell->neighbor(0);
            neighbor->get_dof_indices(local_neighbor_dof_indices);
			local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
			V_neighbor = Cell[local_index].measure();
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                C_R3(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                C_R3(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                C_R3(0,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                C_R3(0,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                C_R3(0,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
            }

            // Row 2 (N neighbour)

			neighbor = cell->neighbor(3);
            neighbor->get_dof_indices(local_neighbor_dof_indices);
			local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
			V_neighbor = Cell[local_index].measure();
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                C_R3(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                C_R3(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                C_R3(1,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                C_R3(1,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                C_R3(1,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
            }
            
            // Row 3 (E neighbour)

			neighbor = cell->neighbor(1);
            neighbor->get_dof_indices(local_neighbor_dof_indices);
			local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
			V_neighbor = Cell[local_index].measure();
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                C_R3(2,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                C_R3(2,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                C_R3(2,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                C_R3(2,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                C_R3(2,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
            }

            // Row 4 (S neighbour) 
            
			neighbor = cell->neighbor(2);
            neighbor->get_dof_indices(local_neighbor_dof_indices);
			local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
			V_neighbor = Cell[local_index].measure();
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                C_R3(3,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                C_R3(3,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                C_R3(3,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                C_R3(3,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                C_R3(3,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
            }

			// Fill A_R3 (Least squares r = 3 matrix)

			ROWS = cell_diagonal_neighbor_index[c].size();
			
			A_R3.reinit(ROWS, 5); A_R3 = 0.0; index = 0; 
			
			// First fill the vertex neighbors 

			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			
	
				local_index = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];
    	        V_neighbor = Cell[local_index].measure();
				
				for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    A_R3(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    A_R3(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    A_R3(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    A_R3(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    A_R3(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
                
                index++; 
			}
			
			CLS_R3[c].initialize(A_R3, C_R3);
            
            // =====================================================================
            // r = 3 stencil 1
            // =====================================================================
            
            if (WW) {

				if (cell_neighbor_index[c][0].size() >= 1) {
					is_admissible_R31[c] = true;
				}
				
				else {
					is_admissible_R31[c] = false;
				}
            }

            else {
				is_admissible_R31[c] = false;
			}
            
            if (is_admissible_R31[c]) {
            
                // Constraint Matrix 
                
                C_R31.reinit(4,5); 
            
                C_R31 = 0.0; 
            
                // Fill C_R31 (C_R31onstraint r = 4 matrix)

                // Row 1 (W neighbour)

				neighbor = cell->neighbor(0);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
    	        V_neighbor = Cell[local_index].measure();
            
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R31(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R31(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R31(0,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R31(0,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R31(0,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }

                // Row 2 (N neighbour)

				neighbor = cell->neighbor(3);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R31(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R31(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R31(1,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R31(1,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R31(1,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }

                // Row 3 (S neighbour) 
                
            
				neighbor = cell->neighbor(2);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R31(2,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R31(2,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R31(2,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R31(2,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R31(2,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
               
				local_index = WW_local_index;
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R31(3,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R31(3,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R31(3,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R31(3,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R31(3,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
                
                // Least Squares Matrix 
                
                index = 0; 

                ROWS = cell_neighbor_index[c][0].size(); 

                A_R31.reinit(ROWS, 5); A_R31 = 0.0; 

				// vertex neighbor of cell at face 0 

				for (unsigned int d = 0; d < cell_neighbor_index[c][0].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_neighbor_index[c][0][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R31(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R31(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R31(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R31(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R31(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
					}
                
					index++; 
				}
								
				CLS_R31[c].initialize(A_R31, C_R31);


            }
            
            // =====================================================================
            // r = 3 stencil 2
            // =====================================================================
            
			if (NN) {
				
				if (cell_neighbor_index[c][3].size() >= 1) {
					is_admissible_R32[c] = true;
				}
				
				else {
					is_admissible_R32[c] = false;
				}
            }

            else {
				is_admissible_R32[c] = false;
			}
            
            if (is_admissible_R32[c]) {
            
                // Constraint Matrix 
                
                C_R32.reinit(4,5); 
            
                C_R32 = 0.0; 
            
                // Fill C_R31 (C_R31onstraint r = 4 matrix)

                // Row 1 (W neighbour)

				neighbor = cell->neighbor(0);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
            
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R32(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R32(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R32(0,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R32(0,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R32(0,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }

                // Row 2 (N neighbour)

				neighbor = cell->neighbor(3);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R32(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R32(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R32(1,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R32(1,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R32(1,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }

                // Row 3 (E neighbour) 
                
				neighbor = cell->neighbor(1);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R32(2,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R32(2,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R32(2,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R32(2,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R32(2,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
                
                // Row 4 (NN neighbour)
                
				local_index = NN_local_index;
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R32(3,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R32(3,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R32(3,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R32(3,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R32(3,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
                
                // Least Squares Matrix 
                
                index = 0;

                ROWS = cell_neighbor_index[c][3].size(); 

                A_R32.reinit(ROWS, 5); A_R32 = 0.0; 

				// vertex neighbor of cell at face 3
			
				for (unsigned int d = 0; d < cell_neighbor_index[c][3].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_neighbor_index[c][3][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R32(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R32(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R32(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R32(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R32(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
					}
                
					index++; 
				}
				
				CLS_R32[c].initialize(A_R32, C_R32);
            }
            
			// =====================================================================
            // r = 3 stencil 3
            // =====================================================================
            
			if (SS) {
				
				if (cell_neighbor_index[c][2].size() >= 1) {
					is_admissible_R33[c] = true;
				}
				
				else {
					is_admissible_R33[c] = false;
				}
            }

            else {
				is_admissible_R33[c] = false;
			}
            
            if (is_admissible_R33[c]) {
            
                // Constraint Matrix 
                
                C_R33.reinit(4,5); 
            
                C_R33 = 0.0; 
		

                // Row 1 (W neighbour)

				neighbor = cell->neighbor(0);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
            
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R33(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R33(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R33(0,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R33(0,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R33(0,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }

                // Row 2 (E neighbour)
				
				neighbor = cell->neighbor(1);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R33(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R33(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R33(1,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R33(1,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R33(1,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }

                // Row 3 (S neighbour) 
                
				neighbor = cell->neighbor(2);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R33(2,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R33(2,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R33(2,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R33(2,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R33(2,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
                
				// Row 4 (SS neighbour) 

				local_index = SS_local_index;
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R33(3,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R33(3,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R33(3,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R33(3,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R33(3,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
                
                // Least Squares Matrix 
                
                index = 0; 

                ROWS = cell_neighbor_index[c][2].size(); 

                A_R33.reinit(ROWS, 5); A_R33 = 0.0; 

				// vertex neighbor of cell at face 2
			
				for (unsigned int d = 0; d < cell_neighbor_index[c][2].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_neighbor_index[c][2][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R33(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R33(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R33(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R33(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R33(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
					}
                
					index++; 
				}
				
                CLS_R33[c].initialize(A_R33, C_R33);
            }
            
            // =====================================================================
            // r = 3 stencil 4
            // =====================================================================
            
			if (EE) {
				
				if (cell_neighbor_index[c][1].size() >= 1) {
					is_admissible_R34[c] = true;
				}
				
				else {
					is_admissible_R34[c] = false;
				}
            }

            else {
				is_admissible_R34[c] = false;
			}
            
            if (is_admissible_R34[c]) {
            
                // Constraint Matrix 
                
                C_R34.reinit(4,5); 
            
                C_R34 = 0.0; 
				
				// Row 1 (N neighbour)

				neighbor = cell->neighbor(3);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
            
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R34(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R34(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R34(0,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R34(0,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R34(0,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }

                // Row 2 (E neighbour)

				neighbor = cell->neighbor(1);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R34(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R34(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R34(1,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R34(1,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R34(1,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }

                // Row 3 (S neighbour) 
                
				neighbor = cell->neighbor(2);
    	        neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R34(2,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R34(2,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R34(2,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R34(2,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R34(2,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
                
				// Row 4 (EE neighbour) 
                
				local_index = EE_local_index;
           	 	V_neighbor = Cell[local_index].measure();
                
                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                    q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                    C_R34(3,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                    C_R34(3,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    C_R34(3,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                    C_R34(3,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                    C_R34(3,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                }
                
                // Least Squares Matrix 
                
                index = 0; 

                ROWS = cell_neighbor_index[c][1].size(); 

                A_R34.reinit(ROWS, 5); A_R34 = 0.0; 

				// vertex neighbor of cell at face 1
			
				for (unsigned int d = 0; d < cell_neighbor_index[c][1].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_neighbor_index[c][1][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R34(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R34(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R34(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R34(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R34(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
					}
                
					index++; 
				}
				
				CLS_R34[c].initialize(A_R34, C_R34);
/*
				if(c == 241){
					pcout<<"A_R34"<<"\tcell: "<<g_i<<"\t local index: "<<c<<"\tcell center: "<<cell->center()<<std::endl;
					for(unsigned int f = 0; f < 4;++f)
						pcout<<"face: "<<f<<"\tface center: "<<cell->face(f)->center()<<std::endl;

					diagonal_cell = cell_neighbor_iterator[c][1].begin();
					end_cell = cell_neighbor_iterator[c][1].end();

					for (; diagonal_cell != end_cell; ++diagonal_cell) {			
						neighbor = *diagonal_cell;
						pcout<<"diagonal center: "<<neighbor->center()<<std::endl;
					}
					pcout<<"EE center: "<<buffer_EE->center()<<std::endl;

					if(!cell->neighbor(1)->at_boundary() ) pcout<<"actual EE center: "<<cell->neighbor(1)->neighbor(1)->center()<<std::endl;
					
					for(unsigned int i = 0;i < index; ++i){
						for(unsigned int j = 0; j < 5 ; ++j)
							pcout<<A_R34(i,j)<<"\t";
						pcout<<std::endl;
					}
					pcout<<"C_R34"<<std::endl;
					for(unsigned int i = 0;i < 4; ++i){
						for(unsigned int j = 0; j < 5 ; ++j)
							pcout<<C_R34(i,j)<<"\t";
						pcout<<std::endl;
					}
				}
*/
            }
            
			// =====================================================================
            // r = 2 stencil 1 (P1 and P2)
            // =====================================================================
            
            A_R21(0,0) = C_R3(0,0); A_R21(0,1) = C_R3(0,1); // Row 1 (P1 cell) 
            A_R21(1,0) = C_R3(1,0); A_R21(1,1) = C_R3(1,1); // Row 2 (P2 cell)

            LU_R21[c].initialize(A_R21);

            // =====================================================================
            // r = 2 stencil 2 (P2 and P3)
            // =====================================================================

            A_R22(0,0) = C_R3(1,0); A_R22(0,1) = C_R3(1,1); // Row 1 (P2 cell)
            A_R22(1,0) = C_R3(2,0); A_R22(1,1) = C_R3(2,1); // Row 2 (P3 cell)

            LU_R22[c].initialize(A_R22);

            // =====================================================================
            // r = 2 stencil 3 (P3 and P4)
            // =====================================================================

            A_R23(0,0) = C_R3(2,0); A_R23(0,1) = C_R3(2,1); // Row 1 (P3 cell)
            A_R23(1,0) = C_R3(3,0); A_R23(1,1) = C_R3(3,1); // Row 2 (P4 cell)

            LU_R23[c].initialize(A_R23);

            // =====================================================================
            // r = 2 stencil 4 (P4 and P1)
            // =====================================================================

            A_R24(0,0) = C_R3(3,0); A_R24(0,1) = C_R3(3,1); // Row 1 (P4 cell)
            A_R24(1,0) = C_R3(0,0); A_R24(1,1) = C_R3(0,1); // Row 2 (P1 cell)

            LU_R24[c].initialize(A_R24);
		} // End of interior cell loop 
        
        else {

            bool W = false, E = false, N = false, S = false;
			if(cell->face(0)->at_boundary()) { W = true; }
			if(cell->face(1)->at_boundary()) { E = true; }
			if(cell->face(2)->at_boundary()) { S = true; }
			if(cell->face(3)->at_boundary()) { N = true; }
            
            // Mark cells at the corner  

            if ((W && N) || (N && E) || (E && S) || (S && W)) {
               is_corner_cell[c] = true;  
            }

            else {
                is_corner_cell[c] = false;
            }

            QGauss<2-1> face_quadrature_formula(2);
            FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);
//            FEFaceValues<2> fv_face_values (mapping, fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);

         	QGauss<2-1> face_center_quadrature_formula(1);
   	     	FEFaceValues<2> fv_face_center_values (fv, face_center_quadrature_formula, update_quadrature_points | update_normal_vectors);
//   	     	FEFaceValues<2> fv_face_center_values (mapping, fv, face_center_quadrature_formula, update_quadrature_points | update_normal_vectors);

            Tensor<1,2> face_normal_vector1; // Face normal vector
	        Tensor<1,2> face_normal_vector2; // Face normal vector
	        Tensor<1,2> face_center_normal_vector; // Face normal vector

            double nx1, ny1;   // Face normal vectors
	        double nx2, ny2; 
	        double nxc, nyc; 

            double x_g1, y_g1, x_g2, y_g2;  

            if (!(is_corner_cell[c])) {
                
                // =====================================================================
                // r = 4 stencil 
                // =====================================================================

            
                Tensor<1,2> face_normal_vector3; // Face normal vector      

                double nx[3], ny[3];   // Face normal vectors 

                double xg[3], yg[3];

                unsigned int f; 

                // Find the boundary face 

				if(W) { f = 0; }
				if(E) { f = 1; }
				if(S) { f = 2; }
				if(N) { f = 3; }

                
                QGauss<2-1> face_quadrature_formula2(3);
                FEFaceValues<2> fv_face_values2 (fv, face_quadrature_formula2, update_quadrature_points | update_normal_vectors);
//                FEFaceValues<2> fv_face_values2 (mapping, fv, face_quadrature_formula2, update_quadrature_points | update_normal_vectors);

                fv_face_values2.reinit(cell, f);

                face_normal_vector1 = fv_face_values2.normal_vector(0); 
                nx[0] = face_normal_vector1[0]; ny[0] = face_normal_vector1[1];
                
                face_normal_vector2 = fv_face_values2.normal_vector(1);
                nx[1] = face_normal_vector2[0]; ny[1] = face_normal_vector2[1];

                face_normal_vector3 = fv_face_values2.normal_vector(2);
                nx[2] = face_normal_vector3[0]; ny[2] = face_normal_vector3[1];
                
                
                xg[0] = fv_face_values2.quadrature_point(0)(0); 
                yg[0] = fv_face_values2.quadrature_point(0)(1);
                xg[1] = fv_face_values2.quadrature_point(1)(0);
                yg[1] = fv_face_values2.quadrature_point(1)(1);
                xg[2] = fv_face_values2.quadrature_point(2)(0);
                yg[2] = fv_face_values2.quadrature_point(2)(1);

                // Constraint Matrix  

                unsigned int index = 0; unsigned int ROWS = 0;  
                C_R4.reinit(7, 9);
                //C_R4.reinit(9,9); 
                C_R4 = 0.0; 
        
                if (!W) {
                    
                    // (W neighbour)

					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
						C_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						C_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						C_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						C_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						C_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
						C_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0)-x0)/h3 - WENO_poly_consts[c](5));
						C_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1)-y0)/h3 - WENO_poly_consts[c](6));
						C_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1)-y0)/h3 - WENO_poly_consts[c](7));
						C_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1)-y0)/h3 - WENO_poly_consts[c](8));
                    }

                    index++;

                } 

                if (!N) {
                    
                    // (N neighbour)

					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
						C_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						C_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						C_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						C_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						C_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
						C_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
						C_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
						C_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
						C_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
                    }

                    index++;

                }

                if (!E) {
                    
                    // (E neighbour)

					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
						C_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						C_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						C_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						C_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						C_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
						C_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
						C_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
						C_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
						C_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
                    }

                    index++;

                }

                if (!S) {
                    
                    // (S neighbor)

					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
						C_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						C_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						C_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						C_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						C_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
						C_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
						C_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
						C_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
						C_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
                    }
                }

                // Boundary conditions 

                // First Derivative 

                for (unsigned int i = 0; i < 3; i++) {
                    C_R4(i+3,0) = nx[i]; 
                    C_R4(i+3,1) = ny[i];
                    C_R4(i+3,2) = 2.0*(xg[i]-x0)*nx[i]/h;
                    C_R4(i+3,3) = 2.0*(yg[i]-y0)*ny[i]/h; 
                    C_R4(i+3,4) = nx[i]*(yg[i]-y0)/h + ny[i]*(xg[i]-x0)/h; 
                    C_R4(i+3,5) = 3.0*nx[i]*(xg[i]-x0)*(xg[i]-x0)/h2; 
                    C_R4(i+3,6) = 3.0*ny[i]*(yg[i]-y0)*(yg[i]-y0)/h2; 
                    C_R4(i+3,7) = 2.0*nx[i]*(xg[i]-x0)*(yg[i]-y0)/h2 + ny[i]*(xg[i]-x0)*(xg[i]-x0)/h2; 
                    C_R4(i+3,8) = 2.0*ny[i]*(xg[i]-x0)*(yg[i]-y0)/h2 + nx[i]*(yg[i]-y0)*(yg[i]-y0)/h2;
                }

                // Third Derivative  
                
                C_R4(6,0) = 0.0; C_R4(6,1) = 0.0; C_R4(6,2) = 0.0; C_R4(6,3) = 0.0;
                C_R4(6,4) = 0.0; C_R4(6,5) = nx[0]; C_R4(6,6) = ny[0]; C_R4(6,7) = 0.0;
                C_R4(6,8) = 0.0;

//				typedef typename std::set<DoFHandler<2>::active_cell_iterator>::iterator dia_cell_iter;
				
				ROWS = cell_diagonal_neighbor_index[c].size(); 
				
				// neighbor of neighbors 
	
				if (WW) {
					ROWS++; 
				}
	
				if (EE) {
					ROWS++; 
				}
		
				if (SS) {
					ROWS++; 
				}

				if (NN) {
					ROWS++; 
				}
				
				index = 0; 
				
				A_R4.reinit(ROWS, 9); A_R4 = 0.0;
			
				// First fill the vertex neighbors 
				
				for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];
    		        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
    	                q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
    	                A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
    	                A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
    	                A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
    	                A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
    	                A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
    	                A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
    	                A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
    	                A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
    	                A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
    	            }
    	            
    	            index++; 
				}

				
				if (WW) {
					
					local_index = WW_local_index;
		            V_neighbor = Cell[local_index].measure();
					
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
						A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
						A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
						A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
						A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
					}
					
					index++; 
				}
				
				if (NN) {
					
					local_index = NN_local_index;
		            V_neighbor = Cell[local_index].measure();
					
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
						A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
						A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
						A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
						A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
					}
					
					index++; 
				}
				
				if (EE) {
					
					local_index = EE_local_index;
		            V_neighbor = Cell[local_index].measure();
					
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
						A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
						A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
						A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
						A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
					}
					
					index++; 
				}
				
				if (SS) {
					
					local_index = SS_local_index;
		            V_neighbor = Cell[local_index].measure();
					
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R4(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R4(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R4(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R4(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R4(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
						A_R4(index,5) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(0) - x0)/h3 - WENO_poly_consts[c](5));
						A_R4(index,6) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](6));
						A_R4(index,7) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](7));
						A_R4(index,8) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)*(q_point(1) - y0)/h3 - WENO_poly_consts[c](8));
					}
					
					index++; 
				}
                
                CLS_R4[c].initialize(A_R4, C_R4);
/*
				if(c == 1){

					std::cout<<"A_R4"<<"\tcell: "<<c<<"\t local index: "<<c<<"\tcell center: "<<cell->center()<<std::endl;

					for(unsigned int f = 0; f < 4;++f)
						pcout<<"face: "<<f<<"\tface center: "<<cell->face(f)->center()<<std::endl;
	
					for(unsigned int i = 0;i < index; ++i){
						for(unsigned int j = 0; j < 9 ; ++j)
							std::cout<<A_R4(i,j)<<"\t";
						std::cout<<std::endl;
					}
					std::cout<<"C_R4"<<std::endl;
					for(unsigned int i = 0;i < 6; ++i){
						for(unsigned int j = 0; j < 9 ; ++j)
							std::cout<<C_R4(i,j)<<"\t";
						std::cout<<std::endl;
					}
				}
*/
            } // End of fourth order transmissive boundaries
		
            // =====================================================================
            // r = 3 centered stencil  
            // =====================================================================
            
            if (!(is_corner_cell[c])) {

				if(cell->face(0)->at_boundary()) { W = true; }
				if(cell->face(1)->at_boundary()) { E = true; }
				if(cell->face(2)->at_boundary()) { S = true; }
				if(cell->face(3)->at_boundary()) { N = true; }
                
                index = 0;  

                // Constraint matrix  
                C_R3.reinit(5, 5);

                C_R3 = 0.0; 

                // Fill C 

                if (!W) {
                    
                    // W neighbour 

					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R3(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R3(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R3(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R3(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R3(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                    
                    index++; 
                }
                
                if (!N) {
                    
                    // (N neighbour)

					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R3(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R3(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R3(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R3(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R3(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                    
                    index++; 
                }
                
                if (!E) {
                    
                    // (E neighbour)

					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R3(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R3(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R3(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R3(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R3(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                    
                    index++; 
                
                }
                
                if (!S) {
                    
                    // (S neighbor)

					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();                    
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R3(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R3(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R3(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R3(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R3(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                    
                    index++; 
                }
                
                // initialize fv_face_values based on the boundary face 
                
                if (W) {
                    fv_face_values.reinit(cell, 0);
                    fv_face_center_values.reinit(cell, 0);

                }
                
                else if (N) {
                    fv_face_values.reinit(cell, 3);
                    fv_face_center_values.reinit(cell, 3);
                }
                
                else if (E) {
                    fv_face_values.reinit(cell, 1);
                    fv_face_center_values.reinit(cell, 1);
                }
                
                else {
                    fv_face_values.reinit(cell, 2);
                    fv_face_center_values.reinit(cell, 2);
                }
                
                face_normal_vector1 = fv_face_values.normal_vector(0); 
                nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                    
                face_normal_vector2 = fv_face_values.normal_vector(1);
                nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];

                face_center_normal_vector = fv_face_center_values.normal_vector(0);
                nxc = face_center_normal_vector[0]; nyc = face_center_normal_vector[1];
                    
                x_g1 = fv_face_values.quadrature_point(0)(0); 
                y_g1 = fv_face_values.quadrature_point(0)(1);
                x_g2 = fv_face_values.quadrature_point(1)(0);
                y_g2 = fv_face_values.quadrature_point(1)(1);
                
                // Boundary Condition on Gauss Point 1 

                C_R3(index,0) = nx1; 
                C_R3(index,1) = ny1; 
                C_R3(index,2) = 2.0*(x_g1-x0)*nx1/h; 
                C_R3(index,3) = 2.0*(y_g1-y0)*ny1/h;
                C_R3(index,4) = (nx1*(y_g1-y0) + ny1*(x_g1-x0))/h;

                // Boundary Condition on Gauss Point 2 

                C_R3(index+1,0) = nx2; 
                C_R3(index+1,1) = ny2; 
                C_R3(index+1,2) = 2.0*(x_g2-x0)*nx2/h; 
                C_R3(index+1,3) = 2.0*(y_g2-y0)*ny2/h;
                C_R3(index+1,4) = (nx2*(y_g2-y0) + ny2*(x_g2-x0))/h; 
                
                CLS_R3[c].initialize(C_R3); 
/*
				if(c == 1){

					std::cout<<"C_R3"<<"\tcell: "<<c<<"\t local index: "<<c<<"\tcell center: "<<cell->center()<<std::endl;

					std::cout<<"C_R3"<<std::endl;
					for(unsigned int i = 0;i < index+1 ; ++i){
						for(unsigned int j = 0; j < 5 ; ++j)
							std::cout<<C_R3(i,j)<<"\t";
						std::cout<<std::endl;
					}
				}
*/				
				// =====================================================================
                // r = 2 stencil 1 
                // =====================================================================
                
                if (W) { // WEST 
                    // P2 cell 

                    A_R21 = 0.0; 
                    
					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R21(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R21(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    A_R21(1,0) = nxc;       A_R21(1,1) = nyc;       // Transmissive boundary
                }
                
                else if (N) { // NORTH 
                    // P1 cell 

                    A_R21 = 0.0; 
                    
					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R21(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R21(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    A_R21(1,0) = nxc;       A_R21(1,1) = nyc;       // Transmissive boundary
                }
                
                else { // WEST and NORTH
                    // P1 cell 

                    A_R21 = 0.0; 
                    
					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R21(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R21(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    // P2 cell 
                    
					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R21(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R21(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                }
                
                LU_R21[c].initialize(A_R21); 
                
                // =====================================================================
                // r = 2 stencil 2
                // =====================================================================
                
                if (E) { // EAST 
                    // P2 cell 

                    A_R22 = 0.0; 
                    
					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R22(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R22(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    A_R22(1,0) = nxc;       A_R22(1,1) = nyc;       // Transmissive boundary
                }
                
                else if (N) { // NORTH 
                    // P3 cell 

                    A_R22 = 0.0; 
                    
					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R22(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R22(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    A_R22(1,0) = nxc;       A_R22(1,1) = nyc;       // Transmissive boundary
                }
                
                else { // NORTH and EAST
                    // P2 cell 

                    A_R22 = 0.0; 
                    
					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R22(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R22(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    // P3 cell 
                    
					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R22(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R22(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                }
                
                LU_R22[c].initialize(A_R22);
                
                // =====================================================================
                // r = 2 stencil 3
                // =====================================================================
                
                if (E) { // EAST 
                    // P4 cell 

                    A_R23 = 0.0; 
                    
					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R23(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R23(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    A_R23(1,0) = nxc;       A_R23(1,1) = nyc;       // Transmissive boundary
                }
                
                else if (S) { // SOUTH 
                    // P3 cell 

                    A_R23 = 0.0; 
                    
					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R23(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R23(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    A_R23(1,0) = nxc;       A_R23(1,1) = nyc;       // Transmissive boundary
                }
                
                else { // EAST and SOUTH 
                    
                    // P3 cell 

                    A_R23 = 0.0; 
                    
					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R23(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R23(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    // P4 cell
                    
					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R23(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R23(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                }
                
                LU_R23[c].initialize(A_R23); 
                
                // =====================================================================
                // r = 2 stencil 4
                // =====================================================================
                
                if (W) { // WEST
                    // P4 cell 

                    A_R24 = 0.0; 
                    
					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R24(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R24(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    A_R24(1,0) = nxc;       A_R24(1,1) = nyc;       // Transmissive boundary
                }
                
                else if (S) { // SOUTH 
                    // P1 cell 

                    A_R24 = 0.0; 
                    
					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R24(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R24(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    A_R24(1,0) = nxc;       A_R24(1,1) = nyc;       // Transmissive boundary
                }
                
                else { // SOUTH and WEST 
                    
                    // P4 cell 

                    A_R24 = 0.0; 
                    
					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R24(0,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R24(0,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                    // P1 cell
                    
					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        A_R24(1,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        A_R24(1,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                    }
                    
                }
                
                LU_R24[c].initialize(A_R24); 
				
            }
            
            // =====================================================================
            // r = 3 stencil 1
            // =====================================================================
            
            if (!W && !(is_corner_cell[c])) {
				
				if (WW) {
					
					if (cell_neighbor_index[c][0].size() >= 1) {
						is_admissible_R31[c] = true;
					}
				
					else {
						is_admissible_R31[c] = false;
					}
					
				}
				
				else {
					is_admissible_R31[c] = false;
				}
            }
            
            else {
                is_admissible_R31[c] = false;
            }
            
            if (is_admissible_R31[c]) {
				
				ROWS = 0; index = 0; 
                
                if (N || S) {
                    C_R31.reinit(5,5); 
                }
        
                else {
                    C_R31.reinit(4,5);
                }
                
                C_R31 = 0.0; 
                
                // Row 0
                
                if (!W) {
            
                    // W neighbour 

					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R31(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R31(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R31(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R31(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R31(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
            
                    index++; 
                }
                
                if (!N) {
                
                    // N neighbour 
                
					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R31(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R31(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R31(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R31(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R31(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (!S) {
                
                    // P4 cell 
                
					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R31(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R31(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R31(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R31(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R31(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (WW) {
                
                    // WW neighbour 
                
					local_index = WW_local_index;
	           	 	V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        C_R31(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R31(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R31(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R31(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R31(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (cell->neighbor_index(2) == -1 || cell->neighbor_index(3) == -1) {
                
                    if (N) {
                        fv_face_values.reinit(cell, 3);
                    }
            
                    else {
                        fv_face_values.reinit(cell, 2);
                    }
                    
                    face_normal_vector1 = fv_face_values.normal_vector(0); 
                    nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                        
                    face_normal_vector2 = fv_face_values.normal_vector(1);
                    nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
                        
                    x_g1 = fv_face_values.quadrature_point(0)(0); 
                    y_g1 = fv_face_values.quadrature_point(0)(1);
                    x_g2 = fv_face_values.quadrature_point(1)(0);
                    y_g2 = fv_face_values.quadrature_point(1)(1);
                    
                    // Boundary Condition on Gauss Point 1 

                    C_R31(index,0) = nx1; 
                    C_R31(index,1) = ny1; 
                    C_R31(index,2) = 2.0*(x_g1-x0)*nx1/h; 
                    C_R31(index,3) = 2.0*(y_g1-y0)*ny1/h;
                    C_R31(index,4) = (nx1*(y_g1-y0) + ny1*(x_g1-x0))/h;

                    // Boundary Condition on Gauss Point 2 

                    C_R31(index+1,0) = nx2; 
                    C_R31(index+1,1) = ny2; 
                    C_R31(index+1,2) = 2.0*(x_g2-x0)*nx2/h; 
                    C_R31(index+1,3) = 2.0*(y_g2-y0)*ny2/h;
                    C_R31(index+1,4) = (nx2*(y_g2-y0) + ny2*(x_g2-x0))/h; 
                }
              
                // Least squares matrix 
                
                index = 0;
                
                ROWS = cell_neighbor_index[c][0].size(); 

                A_R31.reinit(ROWS, 5); A_R31 = 0.0; 
				
				// vertex neighbor of cell at face 0 
			
				for (unsigned int d = 0; d < cell_neighbor_index[c][0].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_neighbor_index[c][0][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R31(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R31(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R31(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R31(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R31(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
					}
                
					index++; 
				}
                
                CLS_R31[c].initialize(A_R31, C_R31);

            } // End of second third order stencil
            
            // =====================================================================
            // r = 3 stencil 2
            // =====================================================================
            
            if (!N && !(is_corner_cell[c])) {
            
				if (NN) {
					
					if (cell_neighbor_index[c][3].size() >= 1) {
						is_admissible_R32[c] = true;
					}
				
					else {
						is_admissible_R32[c] = false;
					}
					
				}
				
				else {
					is_admissible_R32[c] = false;
				}

            }
            
            else {
                is_admissible_R32[c] = false;
            }
            
            if (is_admissible_R32[c]) {
            
                ROWS = 0; index = 0; 

                if (W || E) {
                    C_R32.reinit(5,5); 
                }
        
                else {
                    C_R32.reinit(4,5);
                }
                
                C_R32 = 0.0; 
                
                // Row 0
                
                if (!W) {
            
                    // W cell 

					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        C_R32(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R32(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R32(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R32(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R32(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
            
                    index++; 
                }
                
                if (!N) {
                
                    // P2 cell 
                
					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        C_R32(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R32(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R32(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R32(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R32(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (!E) {
                
                    // P3 cell 
                
					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        C_R32(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R32(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R32(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R32(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R32(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (NN) {
                
                    // S3 cell 
                
					local_index = NN_local_index;
	           	 	V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        C_R32(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R32(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R32(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R32(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R32(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (cell->neighbor_index(0) == -1 || cell->neighbor_index(1) == -1) {
                
                    if (W) {
                        fv_face_values.reinit(cell, 0);
                    }
            
                    else {
                        fv_face_values.reinit(cell, 1);
                    }
                    
                    face_normal_vector1 = fv_face_values.normal_vector(0); 
                    nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                        
                    face_normal_vector2 = fv_face_values.normal_vector(1);
                    nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
                        
                    x_g1 = fv_face_values.quadrature_point(0)(0); 
                    y_g1 = fv_face_values.quadrature_point(0)(1);
                    x_g2 = fv_face_values.quadrature_point(1)(0);
                    y_g2 = fv_face_values.quadrature_point(1)(1);
                    
                    // Boundary Condition on Gauss Point 1 

                    C_R32(index,0) = nx1; 
                    C_R32(index,1) = ny1; 
                    C_R32(index,2) = 2.0*(x_g1-x0)*nx1/h; 
                    C_R32(index,3) = 2.0*(y_g1-y0)*ny1/h;
                    C_R32(index,4) = (nx1*(y_g1-y0) + ny1*(x_g1-x0))/h;

                    // Boundary Condition on Gauss Point 2 

                    C_R32(index+1,0) = nx2; 
                    C_R32(index+1,1) = ny2; 
                    C_R32(index+1,2) = 2.0*(x_g2-x0)*nx2/h; 
                    C_R32(index+1,3) = 2.0*(y_g2-y0)*ny2/h;
                    C_R32(index+1,4) = (nx2*(y_g2-y0) + ny2*(x_g2-x0))/h; 
                }
                
                // Least squares matrix 
                
                index = 0; 

                ROWS = cell_neighbor_index[c][3].size(); 
                A_R32.reinit(ROWS, 5); A_R32 = 0.0; 

				// vertex neighbor of cell at face 3
			
				for (unsigned int d = 0; d < cell_neighbor_index[c][3].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_neighbor_index[c][3][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R32(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R32(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R32(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R32(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R32(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
					}
                
					index++; 
				}				
                
                CLS_R32[c].initialize(A_R32, C_R32);
            } // End of fourth third order stencil transmissive boundaries 
            
            
            // =====================================================================
            // r = 3 stencil 3
            // =====================================================================
            
            if (!S && !(is_corner_cell[c])) {
            
				if (SS) {
					
					if (cell_neighbor_index[c][2].size() >= 1) {
						is_admissible_R33[c] = true;
					}
				
					else {
						is_admissible_R33[c] = false;
					}
					
				}
				
				else {
					is_admissible_R33[c] = false;
				}

            }
            
            else {
                is_admissible_R33[c] = false;
            }
            
            if (is_admissible_R33[c]) {
				
				ROWS = 0; index = 0; 
                
                if (W || E) {
                    C_R33.reinit(5,5); 
                }
        
                else {
                    C_R33.reinit(4,5);
                }
                
                C_R33 = 0.0; 
                
                // Row 0
                
                if (!W) {
            
                    // W cell 

					neighbor = cell->neighbor(0);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R33(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R33(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R33(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R33(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R33(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
            
                    index++; 
                }
                
                if (!E) {
                
                    // E cell 
                
					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R33(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R33(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R33(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R33(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R33(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (!S) {
                
                    // S cell 
                
					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R33(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R33(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R33(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R33(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R33(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (SS) {
                
                    // SS cell 
                
					local_index = SS_local_index;
	           	 	V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R33(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R33(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R33(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R33(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R33(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (cell->neighbor_index(0) == -1 || cell->neighbor_index(1) == -1) {
                
                    if (W) {
                        fv_face_values.reinit(cell, 0);
                    }
            
                    else {
                        fv_face_values.reinit(cell, 1);
                    }
                    
                    face_normal_vector1 = fv_face_values.normal_vector(0); 
                    nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                        
                    face_normal_vector2 = fv_face_values.normal_vector(1);
                    nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
                        
                    x_g1 = fv_face_values.quadrature_point(0)(0); 
                    y_g1 = fv_face_values.quadrature_point(0)(1);
                    x_g2 = fv_face_values.quadrature_point(1)(0);
                    y_g2 = fv_face_values.quadrature_point(1)(1);
                    
                    // Boundary Condition on Gauss Point 1 

                    C_R33(index,0) = nx1; 
                    C_R33(index,1) = ny1; 
                    C_R33(index,2) = 2.0*(x_g1-x0)*nx1/h; 
                    C_R33(index,3) = 2.0*(y_g1-y0)*ny1/h;
                    C_R33(index,4) = (nx1*(y_g1-y0) + ny1*(x_g1-x0))/h;

                    // Boundary Condition on Gauss Point 2 

                    C_R33(index+1,0) = nx2; 
                    C_R33(index+1,1) = ny2; 
                    C_R33(index+1,2) = 2.0*(x_g2-x0)*nx2/h; 
                    C_R33(index+1,3) = 2.0*(y_g2-y0)*ny2/h;
                    C_R33(index+1,4) = (nx2*(y_g2-y0) + ny2*(x_g2-x0))/h; 
                }
                
                // Least squares matrix 
                
				index = 0; 

                ROWS = cell_neighbor_index[c][2].size(); 

				A_R33.reinit(ROWS, 5); A_R33 = 0.0; 

				// vertex neighbor of cell at face 2
			
				for (unsigned int d = 0; d < cell_neighbor_index[c][2].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_neighbor_index[c][2][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R33(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R33(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R33(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R33(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R33(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
					}
                
					index++; 
				}
                
                CLS_R33[c].initialize(A_R33, C_R33);
            } // End of third third order stencil transmissive boundaries 
            
             
            // =====================================================================
            // r = 3 stencil 4
            // =====================================================================
            
            if (!E && !(is_corner_cell[c])) {
            
				if (EE) {
					
					if (cell_neighbor_index[c][1].size() >= 1) {
						is_admissible_R34[c] = true;
					}
				
					else {
						is_admissible_R34[c] = false;
					}
					
				}
				
				else {
					is_admissible_R34[c] = false;
				}

            }
            
            else {
                is_admissible_R34[c] = false;
            }
            
            if (is_admissible_R34[c]) {
				ROWS = 0; index = 0; 
                
                if (N || S) {
                    C_R34.reinit(5,5); 
                }
        
                else {
                    C_R34.reinit(4,5);
                }
                
                C_R34 = 0.0; 
                
                // Row 0
                
                if (!N) {
            
                    // N cell 

					neighbor = cell->neighbor(3);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
                
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R34(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R34(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R34(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R34(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R34(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
            
                    index++; 
                }
                
                if (!E) {
                
                    // P3 cell 
                
					neighbor = cell->neighbor(1);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
                        C_R34(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R34(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R34(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R34(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R34(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (!S) {
                
                    // S cell 
                
					neighbor = cell->neighbor(2);
    		        neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
           	 		V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
				j_w = Cell[local_index].jxw(i);
						C_R34(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R34(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R34(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R34(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R34(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (EE) {
                
                    // EE cell 
                
					local_index = EE_local_index;
	           	 	V_neighbor = Cell[local_index].measure();
            
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
                        C_R34(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
                        C_R34(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
                        C_R34(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
                        C_R34(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
                        C_R34(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
                    }
                
                    index++; 
                }
                
                if (cell->neighbor_index(2) == -1 || cell->neighbor_index(3) == -1) {
                
                    if (N) {
                        fv_face_values.reinit(cell, 3);
                    }
            
                    else {
                        fv_face_values.reinit(cell, 2);
                    }
                    
                    face_normal_vector1 = fv_face_values.normal_vector(0); 
                    nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                        
                    face_normal_vector2 = fv_face_values.normal_vector(1);
                    nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
                        
                    x_g1 = fv_face_values.quadrature_point(0)(0); 
                    y_g1 = fv_face_values.quadrature_point(0)(1);
                    x_g2 = fv_face_values.quadrature_point(1)(0);
                    y_g2 = fv_face_values.quadrature_point(1)(1);
                    
                    // Boundary Condition on Gauss Point 1 

                    C_R34(index,0) = nx1; 
                    C_R34(index,1) = ny1; 
                    C_R34(index,2) = 2.0*(x_g1-x0)*nx1/h; 
                    C_R34(index,3) = 2.0*(y_g1-y0)*ny1/h;
                    C_R34(index,4) = (nx1*(y_g1-y0) + ny1*(x_g1-x0))/h;

                    // Boundary Condition on Gauss Point 2 

                    C_R34(index+1,0) = nx2; 
                    C_R34(index+1,1) = ny2; 
                    C_R34(index+1,2) = 2.0*(x_g2-x0)*nx2/h; 
                    C_R34(index+1,3) = 2.0*(y_g2-y0)*ny2/h;
                    C_R34(index+1,4) = (nx2*(y_g2-y0) + ny2*(x_g2-x0))/h; 
                }
                
                // Least squares matrix 
                
				index = 0; 

                ROWS = cell_neighbor_index[c][1].size(); 

                A_R34.reinit(ROWS, 5); A_R34 = 0.0; 

				// vertex neighbor of cell at face 1
			
				for (unsigned int d = 0; d < cell_neighbor_index[c][1].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_neighbor_index[c][1][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
					for (unsigned int i = 0; i < N_gp*N_gp; i++) {
						q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						A_R34(index,0) += (1./V_neighbor)*j_w*(q_point(0) - WENO_poly_consts[c](0))/h;
						A_R34(index,1) += (1./V_neighbor)*j_w*(q_point(1) - WENO_poly_consts[c](1))/h;
						A_R34(index,2) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(0)-x0)/h2 - WENO_poly_consts[c](2));
						A_R34(index,3) += (1./V_neighbor)*j_w*((q_point(1)-y0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](3));
						A_R34(index,4) += (1./V_neighbor)*j_w*((q_point(0)-x0)*(q_point(1)-y0)/h2 - WENO_poly_consts[c](4));
					}
                
					index++; 
				}
				               
                CLS_R34[c].initialize(A_R34, C_R34);
            } // End of four third order stencil transmissive boundaries 
			
        } // End of boundary cell loop

    } // End of cell loop 
    
    pcout << "Done!" << std::endl;
/*
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "timer.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Time taken to precompute matrices = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}
*/
} // End of function  
