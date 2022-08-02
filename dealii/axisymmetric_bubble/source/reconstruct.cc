#include "../include/Weno432.h"


// Perform the actual reconstruction 

void Weno4_2D::reconstruct() {
	
    Vector<double> W0(5), W1(5), W2(5); 

    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;
    Point<2> face_quadrature_point_1;
    Point<2> face_quadrature_point_2;
/*
	// Gauss Quadrature 
    unsigned int N_gp = 2;  // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);
    Point<2> q_point;
  
    QGauss<2-1> face_quadrature_formula(2);
    FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);
*/
//    Tensor<1,2> face_normal_vector; // Face normal vector

//    double nx, ny;   // Face normal vectors
     
    unsigned int no_stencils = 10;
    unsigned int p = 4; 
    double epsilon = 1.0e-12; 
    double h; // Measure of cell size 
    
    double fourth_order_wt = 0.85;
	double third_order_wt = 0.03; 
	double second_order_wt = 0.0; 
//	double second_order_wt = 0.0125; 
	double sum_gamma; 
	
	Vector<double> gamma(no_stencils);
	
	gamma(0) = fourth_order_wt; 
	gamma(1) = third_order_wt; 
	gamma(2) = third_order_wt;
	gamma(3) = third_order_wt;
	gamma(4) = third_order_wt;
	gamma(5) = third_order_wt;
	gamma(6) = second_order_wt;
	gamma(7) = second_order_wt;
	gamma(8) = second_order_wt;
	gamma(9) = second_order_wt;
    
    // Variables for reconstruction of RHO
    
	double rho0; 
    
	/* Fourth order stencil */ 
	Vector<double> d_rho_4(4); 
	Vector<double> b_rho_4; 
	Vector<double> rho_coeff_4(9); 
    
	/* Third order centered stencil */ 
	Vector<double> d_rho_3(4); 
	Vector<double> b_rho_3; 
	Vector<double> rho_coeff_3(5); 
    
	/* First third order one-sided stencil */ 
	Vector<double> b_rho_31;         
	Vector<double> d_rho_31;
	Vector<double> rho_coeff_31(5);
    
	/* Second third order one-sided stencil */
	Vector<double> b_rho_32;       
	Vector<double> d_rho_32;
	Vector<double> rho_coeff_32(5);
    
    /* Third third order one-sided stencil */
	Vector<double> b_rho_33;          
	Vector<double> d_rho_33;
	Vector<double> rho_coeff_33(5);
    
    /* Fourth third order one-sided stencil */
    Vector<double> b_rho_34; 
    Vector<double> d_rho_34;
    Vector<double> rho_coeff_34(5);  
	
	Vector<double> rho_coeff_21(2); 
	Vector<double> rho_coeff_22(2); 
	Vector<double> rho_coeff_23(2); 
	Vector<double> rho_coeff_24(2);
	Vector<double> b_rho2(2);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_RHO(no_stencils); Vector<double> w_RHO(no_stencils); double sum_RHO;
    
    // Variables for reconstruction of U
    // Variables for reconstruction of U
    double u0; 

	/* Fourth order stencil */ 
	Vector<double> d_u_4(4); 
	Vector<double> b_u_4; 
	Vector<double> u_coeff_4(9);
    
	/* Third order centered stencil */ 
	Vector<double> d_u_3(4); 
	Vector<double> b_u_3; 
	Vector<double> u_coeff_3(5); 

	/* First third order one-sided stencil */ 
	Vector<double> b_u_31;         
	Vector<double> d_u_31;
	Vector<double> u_coeff_31(5);
    
	/* Second third order one-sided stencil */
	Vector<double> b_u_32;       
	Vector<double> d_u_32;
	Vector<double> u_coeff_32(5);
    
	/* Third third order one-sided stencil */
	Vector<double> b_u_33;          
	Vector<double> d_u_33;
	Vector<double> u_coeff_33(5);
    
	/* Fourth third order one-sided stencil */
    Vector<double> b_u_34; 
    Vector<double> d_u_34;
    Vector<double> u_coeff_34(5);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_U(no_stencils); Vector<double> w_U(no_stencils);  double sum_U;
	
	Vector<double> u_coeff_21(2); 
    Vector<double> u_coeff_22(2); 
    Vector<double> u_coeff_23(2); 
    Vector<double> u_coeff_24(2); 
    Vector<double> b_u_2(2);
    
    // Variables for reconstruction of V
    // Variables for reconstruction of V
    double v0; 

	/* Fourth order stencil */ 
	Vector<double> d_v_4(4); 
	Vector<double> b_v_4; 
	Vector<double> v_coeff_4(9);
    
	/* Third order centered stencil */ 
	Vector<double> d_v_3(4); 
	Vector<double> b_v_3; 
	Vector<double> v_coeff_3(5); 

	/* First third order one-sided stencil */ 
	Vector<double> b_v_31;         
	Vector<double> d_v_31;
	Vector<double> v_coeff_31(5);
    
	/* Second third order one-sided stencil */
	Vector<double> b_v_32;       
	Vector<double> d_v_32;
	Vector<double> v_coeff_32(5);
    
	/* Third third order one-sided stencil */
	Vector<double> b_v_33;          
	Vector<double> d_v_33;
	Vector<double> v_coeff_33(5);
    
	/* Fourth third order one-sided stencil */
    Vector<double> b_v_34; 
    Vector<double> d_v_34;
    Vector<double> v_coeff_34(5);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_V(no_stencils); Vector<double> w_V(no_stencils);  double sum_V;
	
	Vector<double> v_coeff_21(2); 
    Vector<double> v_coeff_22(2); 
    Vector<double> v_coeff_23(2); 
    Vector<double> v_coeff_24(2); 
    Vector<double> b_v_2(2);
    
    // Variables for reconstruction of E
    // Variables for reconstruction of P
    double p0; 

	/* Fourth order stencil */ 
	Vector<double> d_p_4(4); 
	Vector<double> b_p_4; 
	Vector<double> p_coeff_4(9);
    
	/* Third order centered stencil */ 
	Vector<double> d_p_3(4); 
	Vector<double> b_p_3; 
	Vector<double> p_coeff_3(5); 

	/* First third order one-sided stencil */ 
	Vector<double> b_p_31;         
	Vector<double> d_p_31;
	Vector<double> p_coeff_31(5);
    
	/* Second third order one-sided stencil */
	Vector<double> b_p_32;       
	Vector<double> d_p_32;
	Vector<double> p_coeff_32(5);
    
	/* Third third order one-sided stencil */
	Vector<double> b_p_33;          
	Vector<double> d_p_33;
	Vector<double> p_coeff_33(5);
    
	/* Fourth third order one-sided stencil */
    Vector<double> b_p_34; 
    Vector<double> d_p_34;
    Vector<double> p_coeff_34(5);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_P(no_stencils); Vector<double> w_P(no_stencils);  double sum_P;
	
	Vector<double> p_coeff_21(2); 
    Vector<double> p_coeff_22(2); 
    Vector<double> p_coeff_23(2); 
    Vector<double> p_coeff_24(2); 
    Vector<double> b_p_2(2);


    // Variables for reconstruction of phi
    double phi0; 

	/* Fourth order stencil */ 
	Vector<double> d_phi_4(4); 
	Vector<double> b_phi_4; 
	Vector<double> phi_coeff_4(9);
    
	/* Third order centered stencil */ 
	Vector<double> d_phi_3(4); 
	Vector<double> b_phi_3; 
	Vector<double> phi_coeff_3(5); 

	/* First third order one-sided stencil */ 
	Vector<double> b_phi_31;         
	Vector<double> d_phi_31;
	Vector<double> phi_coeff_31(5);
    
	/* Second third order one-sided stencil */
	Vector<double> b_phi_32;       
	Vector<double> d_phi_32;
	Vector<double> phi_coeff_32(5);
    
	/* Third third order one-sided stencil */
	Vector<double> b_phi_33;          
	Vector<double> d_phi_33;
	Vector<double> phi_coeff_33(5);
    
	/* Fourth third order one-sided stencil */
    Vector<double> b_phi_34; 
    Vector<double> d_phi_34;
    Vector<double> phi_coeff_34(5);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_PHI(no_stencils); Vector<double> w_PHI(no_stencils);  double sum_PHI;
	
	Vector<double> phi_coeff_21(2); 
    Vector<double> phi_coeff_22(2); 
    Vector<double> phi_coeff_23(2); 
    Vector<double> phi_coeff_24(2); 
    Vector<double> b_phi_2(2);
    
    // Iterate over all the cells 
    
    DoFHandler<2>::active_cell_iterator cell, neighbor;

	double gamma_ , gamma_1, gamma_2, p_inf_, p_inf_1, p_inf_2;
    
    unsigned int index, ROWS, g_i; 

	unsigned int WW_index, NN_index, EE_index, SS_index;

	unsigned int neighbor_p1, neighbor_p2, neighbor_p3, neighbor_p4;

//   	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)	
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {
		//std::cout<<"c: "<<c<<std::endl;
		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];

		neighbor_p1 = numbers::invalid_unsigned_int;
		neighbor_p2 = numbers::invalid_unsigned_int;
		neighbor_p3 = numbers::invalid_unsigned_int;
		neighbor_p4 = numbers::invalid_unsigned_int;

		WW_index = numbers::invalid_unsigned_int;
		NN_index = numbers::invalid_unsigned_int;
		EE_index = numbers::invalid_unsigned_int;
		SS_index = numbers::invalid_unsigned_int;

	    bool WW = false, EE = false, NN = false, SS = false;

        rho0	=   Rho(c);
        u0		=   U(c);
        v0 		=   V(c);
        p0     	=   P(c);
        phi0    =   Phi(c);
        h = std::sqrt(Cell[c].measure()); 

		if(!cell->face(0)->at_boundary()) {
	        neighbor = cell->neighbor(0);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p1 = global_to_local_index_map[local_neighbor_dof_indices[0] ];
		}

		if(!cell->face(3)->at_boundary()) {
	        neighbor = cell->neighbor(3);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p2 = global_to_local_index_map[local_neighbor_dof_indices[0] ];
		}

		if(!cell->face(1)->at_boundary()) {
	        neighbor = cell->neighbor(1);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p3 = global_to_local_index_map[local_neighbor_dof_indices[0] ];
		}

		if(!cell->face(2)->at_boundary()) {
	        neighbor = cell->neighbor(2);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p4 = global_to_local_index_map[local_neighbor_dof_indices[0] ];
		}

		if (cell_neighbor_neighbor_index[c][0].size() > 0) {
			WW = true;
			WW_index = global_to_local_index_map[cell_neighbor_neighbor_index[c][0][0] ];
		}

		if (cell_neighbor_neighbor_index[c][1].size() > 0) {
			EE = true;
			EE_index = global_to_local_index_map[cell_neighbor_neighbor_index[c][1][0] ];
		}

		if (cell_neighbor_neighbor_index[c][2].size() > 0) {
			SS = true;
			SS_index = global_to_local_index_map[cell_neighbor_neighbor_index[c][2][0] ];
		}

		if (cell_neighbor_neighbor_index[c][3].size() > 0) {
			NN = true;
			NN_index = global_to_local_index_map[cell_neighbor_neighbor_index[c][3][0] ];
		}
        
        if ( !(cell->at_boundary()) ) {
//			//std::cout<<"interior"<<std::endl;
					//std::cout<<"c in: "<<c<<std::endl;
            d_rho_31.reinit(4);   d_rho_32.reinit(4);   d_rho_33.reinit(4);   d_rho_34.reinit(4);
            d_u_31.reinit(4); d_u_32.reinit(4); d_u_33.reinit(4); d_u_34.reinit(4);
            d_v_31.reinit(4); d_v_32.reinit(4); d_v_33.reinit(4); d_v_34.reinit(4);
            d_p_31.reinit(4);     d_p_32.reinit(4);     d_p_33.reinit(4);     d_p_34.reinit(4);
            d_phi_31.reinit(4);     d_phi_32.reinit(4);     d_phi_33.reinit(4);     d_phi_34.reinit(4);
            
            // =====================================================================
            // r = 4 stencil 
            // =====================================================================
            
			d_rho_4(0)  = (Rho(neighbor_p1) - rho0);   // W neighbor 
			d_rho_4(1)  = (Rho(neighbor_p2) - rho0);   // N neighbor
			d_rho_4(2)  = (Rho(neighbor_p3) - rho0);   // E neighbor
			d_rho_4(3)  = (Rho(neighbor_p4) - rho0);   // S neighbor
			
			d_u_4(0)  = (U(neighbor_p1) - u0);   // W neighbor 
			d_u_4(1)  = (U(neighbor_p2) - u0);    // N neighbor
			d_u_4(2)  = (U(neighbor_p3) - u0);    // E neighbor
			d_u_4(3)  = (U(neighbor_p4) - u0);    // S neighbor
			
			d_v_4(0)  = (V(neighbor_p1) - v0);   // W neighbor 
			d_v_4(1)  = (V(neighbor_p2) - v0);   // N neighbor
			d_v_4(2)  = (V(neighbor_p3) - v0);   // S neighbor
			d_v_4(3)  = (V(neighbor_p4) - v0);   // E neighbor
			
			d_p_4(0)  = (P(neighbor_p1) - p0);   // W neighbor 
			d_p_4(1)  = (P(neighbor_p2) - p0);   // N neighbor
			d_p_4(2)  = (P(neighbor_p3) - p0);   // S neighbor
			d_p_4(3)  = (P(neighbor_p4) - p0);   // E neighbor

			d_phi_4(0)  = (Phi(neighbor_p1) - phi0);   // W neighbor 
			d_phi_4(1)  = (Phi(neighbor_p2) - phi0);   // N neighbor
			d_phi_4(2)  = (Phi(neighbor_p3) - phi0);   // S neighbor
			d_phi_4(3)  = (Phi(neighbor_p4) - phi0);   // E neighbor
            
            index = 0; 
            
            // Least Squares Part  
		
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
            
            b_rho_4.reinit(ROWS); b_u_4.reinit(ROWS); b_v_4.reinit(ROWS); b_p_4.reinit(ROWS); b_phi_4.reinit(ROWS);
			
			// vertex neighbors
			
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];

				b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
				b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
				b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
				b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
				b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
				index++; 
			}
			
			// neighbors of neighbors 
			
			if (WW) {

				local_neighbor_dof_indices[0] = WW_index;

				b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
				b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
				b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
				b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
				b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
				index++; 
			}
			
			if (NN) {

				local_neighbor_dof_indices[0] = NN_index;

				b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
				b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
				b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
				b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
				b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
				index++; 
			}
			
			if (EE) {

				local_neighbor_dof_indices[0] = EE_index;

				b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
				b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
				b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
				b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
				b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
				index++; 
			}
			
			if (SS) {

				local_neighbor_dof_indices[0] = SS_index;

				b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
				b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
				b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
				b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
				b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
				index++; 
			}

			CLS_R4[c].solve(b_rho_4, d_rho_4, rho_coeff_4);  
			CLS_R4[c].solve(b_u_4, d_u_4, u_coeff_4);
			CLS_R4[c].solve(b_v_4, d_v_4, v_coeff_4);
			CLS_R4[c].solve(b_p_4, d_p_4, p_coeff_4);
			CLS_R4[c].solve(b_phi_4, d_phi_4, phi_coeff_4);

			// =====================================================================
            // r = 3 stencil (Centered stencil)
            // =====================================================================
			
			// constraint part (consists of face neighbours)
			
			d_rho_3(0)  = (Rho(neighbor_p1) - rho0);                // W neighbor 
			d_rho_3(1)  = (Rho(neighbor_p2) - rho0);                // N neighbor
			d_rho_3(2)  = (Rho(neighbor_p3) - rho0);                // E neighbor
			d_rho_3(3)  = (Rho(neighbor_p4) - rho0);                // S neighbor
			
			d_u_3(0)  = (U(neighbor_p1) - u0);          // W neighbor 
			d_u_3(1)  = (U(neighbor_p2) - u0);          // N neighbor
			d_u_3(2)  = (U(neighbor_p3) - u0);          // E neighbor
			d_u_3(3)  = (U(neighbor_p4) - u0);          // S neighbor
			
			d_v_3(0)  = (V(neighbor_p1) - v0);          // W neighbor 
			d_v_3(1)  = (V(neighbor_p2) - v0);          // N neighbor
			d_v_3(2)  = (V(neighbor_p3) - v0);          // E neighbor
			d_v_3(3)  = (V(neighbor_p4) - v0);          // S neighbor
			
			d_p_3(0)  = (P(neighbor_p1) - p0);                      // W neighbor 
			d_p_3(1)  = (P(neighbor_p2) - p0);                      // N neighbor
			d_p_3(2)  = (P(neighbor_p3) - p0);                      // E neighbor
			d_p_3(3)  = (P(neighbor_p4) - p0);                      // S neighbor

			d_phi_3(0)  = (Phi(neighbor_p1) - phi0);                      // W neighbor 
			d_phi_3(1)  = (Phi(neighbor_p2) - phi0);                      // N neighbor
			d_phi_3(2)  = (Phi(neighbor_p3) - phi0);                      // E neighbor
			d_phi_3(3)  = (Phi(neighbor_p4) - phi0);                      // S neighbor
					
			ROWS = cell_diagonal_neighbor_index[c].size();
			index = 0; 
            
            b_rho_3.reinit(ROWS);
            b_u_3.reinit(ROWS);
            b_v_3.reinit(ROWS);
            b_p_3.reinit(ROWS);
            b_phi_3.reinit(ROWS);

			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];

				b_rho_3(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
				b_u_3(index) = U(local_neighbor_dof_indices[0] ) - u0;
				b_v_3(index) = V(local_neighbor_dof_indices[0] ) - v0;
				b_p_3(index)     = P(local_neighbor_dof_indices[0] ) - p0;
				b_phi_3(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
				index++; 
			}

			CLS_R3[c].solve(b_rho_3, d_rho_3, rho_coeff_3);
			CLS_R3[c].solve(b_u_3, d_u_3, u_coeff_3); 
			CLS_R3[c].solve(b_v_3, d_v_3, v_coeff_3);
			CLS_R3[c].solve(b_p_3, d_p_3, p_coeff_3);
			CLS_R3[c].solve(b_phi_3, d_phi_3, phi_coeff_3);

			
			// =====================================================================
			// r = 3 stencil 1
			// =====================================================================
            
			if (is_admissible_R31[c]) {
            
				d_rho_31(0)  = Rho(neighbor_p1) - rho0;                // W neighbor 
				d_rho_31(1)  = Rho(neighbor_p2) - rho0;                // N neighbor
				d_rho_31(2)  = Rho(neighbor_p4) - rho0;                // S neighbor
				d_rho_31(3)  = Rho(WW_index) - rho0;               // WW neighbor
                
				d_u_31(0)  = U(neighbor_p1) - u0;          // W neighbor 
				d_u_31(1)  = U(neighbor_p2) - u0;          // N neighbor
				d_u_31(2)  = U(neighbor_p4) - u0;          // S neighbor
				d_u_31(3)  = U(WW_index) - u0;         // WW neighbor
                
				d_v_31(0)  = V(neighbor_p1) - v0;          // W neighbor 
				d_v_31(1)  = V(neighbor_p2) - v0;          // N neighbor
				d_v_31(2)  = V(neighbor_p4) - v0;          // S neighbor
				d_v_31(3)  = V(WW_index) - v0;         // WW neighbor
                
				d_p_31(0)  = P(neighbor_p1) - p0;                      // W neighbor 
				d_p_31(1)  = P(neighbor_p2) - p0;                      // N neighbor
				d_p_31(2)  = P(neighbor_p4) - p0;                      // S neighbor
				d_p_31(3)  = P(WW_index) - p0;                     // WW neighbor

				d_phi_31(0)  = Phi(neighbor_p1) - phi0;                      // W neighbor 
				d_phi_31(1)  = Phi(neighbor_p2) - phi0;                      // N neighbor
				d_phi_31(2)  = Phi(neighbor_p4) - phi0;                      // S neighbor
				d_phi_31(3)  = Phi(WW_index) - phi0;                     // WW neighbor
				
                ROWS = cell_neighbor_index[c][0].size(); 
				
				b_rho_31.reinit(ROWS); b_u_31.reinit(ROWS); b_v_31.reinit(ROWS); b_p_31.reinit(ROWS);
				b_phi_31.reinit(ROWS);

                index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][0].size(); ++d) {
	
					local_neighbor_dof_indices[0] = global_to_local_index_map[cell_neighbor_index[c][0][d] ];

					b_rho_31(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_31(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_31(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_31(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_31(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
                CLS_R31[c].solve(b_rho_31, d_rho_31, rho_coeff_31);
                CLS_R31[c].solve(b_u_31, d_u_31, u_coeff_31);
                CLS_R31[c].solve(b_v_31, d_v_31, v_coeff_31);
                CLS_R31[c].solve(b_p_31, d_p_31, p_coeff_31); 
                CLS_R31[c].solve(b_phi_31, d_phi_31, phi_coeff_31); 
			}
            
            // If stencil is not available fallback to third order centered stencil 
            
            else {
				rho_coeff_31 = rho_coeff_3; 
				u_coeff_31 = u_coeff_3;
				v_coeff_31 = v_coeff_3;
				p_coeff_31 = p_coeff_3;
				phi_coeff_31 = phi_coeff_3;
			}
            
            // =====================================================================
            // r = 3 stencil 2
            // =====================================================================
            
			if (is_admissible_R32[c]) {
            
				d_rho_32(0)  = Rho(neighbor_p1) - rho0;                // W neighbor 
				d_rho_32(1)  = Rho(neighbor_p2) - rho0;                // N neighbor
				d_rho_32(2)  = Rho(neighbor_p3) - rho0;                // E neighbor
				d_rho_32(3)  = Rho(NN_index) - rho0;               // NN neighbor
                
				d_u_32(0)  = U(neighbor_p1) - u0;          // W neighbor 
				d_u_32(1)  = U(neighbor_p2) - u0;          // N neighbor
				d_u_32(2)  = U(neighbor_p3) - u0;          // E neighbor
				d_u_32(3)  = U(NN_index) - u0;           // NN neighbor
                
				d_v_32(0)  = V(neighbor_p1) - v0;          // W neighbor 
				d_v_32(1)  = V(neighbor_p2) - v0;          // N neighbor
				d_v_32(2)  = V(neighbor_p3) - v0;          // E neighbor
				d_v_32(3)  = V(NN_index) - v0;           // NN neighbor
                
				d_p_32(0)  = P(neighbor_p1) - p0;                      // W neighbor 
				d_p_32(1)  = P(neighbor_p2) - p0;                      // N neighbor
				d_p_32(2)  = P(neighbor_p3) - p0;                      // E neighbor
				d_p_32(3)  = P(NN_index) - p0;                     // NN neighbor

				d_phi_32(0)  = Phi(neighbor_p1) - phi0;                      // W neighbor 
				d_phi_32(1)  = Phi(neighbor_p2) - phi0;                      // N neighbor
				d_phi_32(2)  = Phi(neighbor_p3) - phi0;                      // E neighbor
				d_phi_32(3)  = Phi(NN_index) - phi0;                     // NN neighbor

                ROWS = cell_neighbor_index[c][3].size(); 

				b_rho_32.reinit(ROWS); b_u_32.reinit(ROWS); b_v_32.reinit(ROWS); b_p_32.reinit(ROWS);
				b_phi_32.reinit(ROWS);

				index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][3].size(); ++d) {
	
					local_neighbor_dof_indices[0] = global_to_local_index_map[cell_neighbor_index[c][3][d] ];

					b_rho_32(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_32(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_32(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_32(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_32(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
                CLS_R32[c].solve(b_rho_32, d_rho_32, rho_coeff_32);
                CLS_R32[c].solve(b_u_32, d_u_32, u_coeff_32);
                CLS_R32[c].solve(b_v_32, d_v_32, v_coeff_32);
                CLS_R32[c].solve(b_p_32, d_p_32, p_coeff_32); 
                CLS_R32[c].solve(b_phi_32, d_phi_32, phi_coeff_32); 
			}
            // If stencil is not available fallback to third order centered stencil
            
            else {
				rho_coeff_32 = rho_coeff_3; 
				u_coeff_32 = u_coeff_3;
				v_coeff_32 = v_coeff_3;
				p_coeff_32 = p_coeff_3;
				phi_coeff_32 = phi_coeff_3;
            }  
            
            // =====================================================================
            // r = 3 stencil 3
            // =====================================================================
            
			if (is_admissible_R33[c]) {
            
				d_rho_33(0)  = Rho(neighbor_p1) - rho0;                // W neighbor 
				d_rho_33(1)  = Rho(neighbor_p3) - rho0;                // E neighbor
				d_rho_33(2)  = Rho(neighbor_p4) - rho0;                // S neighbor
				d_rho_33(3)  = Rho(SS_index) - rho0;               // SS neighbor
                
				d_u_33(0)  = U(neighbor_p1) - u0;          // W neighbor 
				d_u_33(1)  = U(neighbor_p3) - u0;          // E neighbor
				d_u_33(2)  = U(neighbor_p4) - u0;          // S neighbor
				d_u_33(3)  = U(SS_index) - u0;         // SS neighbor
                
				d_v_33(0)  = V(neighbor_p1) - v0;          // W neighbor 
				d_v_33(1)  = V(neighbor_p3) - v0;          // E neighbor
				d_v_33(2)  = V(neighbor_p4) - v0;          // S neighbor
				d_v_33(3)  = V(SS_index) - v0;         // SS neighbor
                
				d_p_33(0)  = P(neighbor_p1) - p0;                      // W neighbor 
				d_p_33(1)  = P(neighbor_p3) - p0;                      // E neighbor
				d_p_33(2)  = P(neighbor_p4) - p0;                      // S neighbor
				d_p_33(3)  = P(SS_index) - p0;                     // SS neighbor

				d_phi_33(0)  = Phi(neighbor_p1) - phi0;                      // W neighbor 
				d_phi_33(1)  = Phi(neighbor_p3) - phi0;                      // E neighbor
				d_phi_33(2)  = Phi(neighbor_p4) - phi0;                      // S neighbor
				d_phi_33(3)  = Phi(SS_index) - phi0;                     // SS neighbor

                ROWS = cell_neighbor_index[c][2].size(); 

				b_rho_33.reinit(ROWS); b_u_33.reinit(ROWS); b_v_33.reinit(ROWS); b_p_33.reinit(ROWS);
				b_phi_33.reinit(ROWS);

                index = 0; 			

				for (unsigned int d = 0; d < cell_neighbor_index[c][2].size(); ++d) {
	
					local_neighbor_dof_indices[0] = global_to_local_index_map[cell_neighbor_index[c][2][d] ];

					b_rho_33(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_33(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_33(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_33(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_33(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
                CLS_R33[c].solve(b_rho_33,   d_rho_33,   rho_coeff_33);
                CLS_R33[c].solve(b_u_33, d_u_33, u_coeff_33);
                CLS_R33[c].solve(b_v_33, d_v_33, v_coeff_33);
                CLS_R33[c].solve(b_p_33,     d_p_33,     p_coeff_33); 
                CLS_R33[c].solve(b_phi_33,     d_phi_33,     phi_coeff_33); 
			}
            // If stencil is not available fallback to third order centered stencil
            
            else {
				rho_coeff_33 = rho_coeff_3; 
				u_coeff_33 = u_coeff_3;
				v_coeff_33 = v_coeff_3;
				p_coeff_33 = p_coeff_3;
				phi_coeff_33 = phi_coeff_3;
			} 
            
			// =====================================================================
            // r = 3 stencil 4
            // =====================================================================
            
			if (is_admissible_R34[c]) {
            
				d_rho_34(0)  = Rho(neighbor_p2) - rho0;                // N neighbor 
				d_rho_34(1)  = Rho(neighbor_p3) - rho0;                // E neighbor
				d_rho_34(2)  = Rho(neighbor_p4) - rho0;                // S neighbor
				d_rho_34(3)  = Rho(EE_index) - rho0;                // S neighbor
                
				d_u_34(0)  = U(neighbor_p2) - u0;          // N neighbor 
				d_u_34(1)  = U(neighbor_p3) - u0;          // E neighbor
				d_u_34(2)  = U(neighbor_p4) - u0;          // S neighbor
				d_u_34(3)  = U(EE_index) - u0;          // S neighbor
                
				d_v_34(0)  = V(neighbor_p2) - v0;          // N neighbor 
				d_v_34(1)  = V(neighbor_p3) - v0;          // E neighbor
				d_v_34(2)  = V(neighbor_p4) - v0;          // S neighbor
				d_v_34(3)  = V(EE_index) - v0;          // S neighbor
                
				d_p_34(0)  = P(neighbor_p2) - p0;                      // N neighbor 
				d_p_34(1)  = P(neighbor_p3) - p0;                      // E neighbor
				d_p_34(2)  = P(neighbor_p4) - p0;                      // S neighbor
				d_p_34(3)  = P(EE_index) - p0;                      // S neighbor

				d_phi_34(0)  = Phi(neighbor_p2) - phi0;                      // N neighbor 
				d_phi_34(1)  = Phi(neighbor_p3) - phi0;                      // E neighbor
				d_phi_34(2)  = Phi(neighbor_p4) - phi0;                      // S neighbor
				d_phi_34(3)  = Phi(EE_index) - phi0;                      // S neighbor

                ROWS = cell_neighbor_index[c][1].size(); 

                b_rho_34.reinit(ROWS); b_u_34.reinit(ROWS); b_v_34.reinit(ROWS); b_p_34.reinit(ROWS);
				b_phi_34.reinit(ROWS);

                index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][1].size(); ++d) {
	
					local_neighbor_dof_indices[0] = global_to_local_index_map[cell_neighbor_index[c][1][d] ];

					b_rho_34(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_34(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_34(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_34(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_34(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
				CLS_R34[c].solve(b_rho_34, d_rho_34, rho_coeff_34);
				CLS_R34[c].solve(b_u_34, d_u_34, u_coeff_34);
				CLS_R34[c].solve(b_v_34, d_v_34, v_coeff_34);
				CLS_R34[c].solve(b_p_34, d_p_34, p_coeff_34); 
				CLS_R34[c].solve(b_phi_34, d_phi_34, phi_coeff_34); 
            }
            
            // If stencil is not available fallback to third order centered stencil
            
            else {
				rho_coeff_34 = rho_coeff_3; 
				u_coeff_34 = u_coeff_3;
				v_coeff_34 = v_coeff_3;
				p_coeff_34 = p_coeff_3;
				phi_coeff_34 = phi_coeff_3;
            } 
            
            
			// =====================================================================
			// r = 2 stencil 1
			// ===================================================================== 

			b_rho2(0)  = (Rho(neighbor_p1) - rho0);          // W neighbor 
			b_rho2(1)  = (Rho(neighbor_p2) - rho0);          // N neighbor
			LU_R21[c].solve(b_rho2, rho_coeff_21);

			b_u_2(0)  = (U(neighbor_p1) - u0);   // W neighbor 
			b_u_2(1)  = (U(neighbor_p2) - u0);   // N neighbor
			LU_R21[c].solve(b_u_2, u_coeff_21);

			b_v_2(0)  = (V(neighbor_p1) - v0);   // W neighbor 
			b_v_2(1)  = (V(neighbor_p2) - v0);   // N neighbor
			LU_R21[c].solve(b_v_2, v_coeff_21);

			b_p_2(0)  = (P(neighbor_p1) - p0);               // W neighbor 
			b_p_2(1)  = (P(neighbor_p2) - p0);               // N neighbor
			LU_R21[c].solve(b_p_2, p_coeff_21);

			b_phi_2(0)  = (Phi(neighbor_p1) - phi0);               // W neighbor 
			b_phi_2(1)  = (Phi(neighbor_p2) - phi0);               // N neighbor
			LU_R21[c].solve(b_phi_2, phi_coeff_21);
            
            // =====================================================================
            // r = 2 stencil 2 
            // ===================================================================== 

			b_rho2(0) = (Rho(neighbor_p2) - rho0);          // N neighbor
			b_rho2(1) = (Rho(neighbor_p3) - rho0);          // E neighbor
			LU_R22[c].solve(b_rho2, rho_coeff_22);

			b_u_2(0) = (U(neighbor_p2) - u0);   // N neighbor
			b_u_2(1) = (U(neighbor_p3) - u0);   // E neighbor
			LU_R22[c].solve(b_u_2, u_coeff_22);

			b_v_2(0) = (V(neighbor_p2) - v0);   // N neighbor
			b_v_2(1) = (V(neighbor_p3) - v0);   // E neighbor
			LU_R22[c].solve(b_v_2, v_coeff_22);

			b_p_2(0) = (P(neighbor_p2) - p0);               // N neighbor
			b_p_2(1) = (P(neighbor_p3) - p0);               // E neighbor
			LU_R22[c].solve(b_p_2, p_coeff_22);

			b_phi_2(0) = (Phi(neighbor_p2) - phi0);               // N neighbor
			b_phi_2(1) = (Phi(neighbor_p3) - phi0);               // E neighbor
			LU_R22[c].solve(b_phi_2, phi_coeff_22);
		
			// =====================================================================
			// r = 2 stencil 3
			// =====================================================================

			b_rho2(0) = (Rho(neighbor_p3) - rho0);          // E neighbor
			b_rho2(1) = (Rho(neighbor_p4) - rho0);          // S neighbor                                          
			LU_R23[c].solve(b_rho2, rho_coeff_23); 

			b_u_2(0) = (U(neighbor_p3) - u0);   // E neighbor
			b_u_2(1) = (U(neighbor_p4) - u0);   // S neighbor                                   
			LU_R23[c].solve(b_u_2, u_coeff_23); 

			b_v_2(0) = (V(neighbor_p3) - v0);   // E neighbor
			b_v_2(1) = (V(neighbor_p4) - v0);   // S neighbor                                   
			LU_R23[c].solve(b_v_2, v_coeff_23);

			b_p_2(0) = (P(neighbor_p3) - p0);               // E neighbor
			b_p_2(1) = (P(neighbor_p4) - p0);               // S neighbor                                              
			LU_R23[c].solve(b_p_2, p_coeff_23);

			b_phi_2(0) = (Phi(neighbor_p3) - phi0);               // E neighbor
			b_phi_2(1) = (Phi(neighbor_p4) - phi0);               // S neighbor                                              
			LU_R23[c].solve(b_phi_2, phi_coeff_23);  
			
			// =====================================================================
			// r = 2 stencil 4
			// =====================================================================

			b_rho2(0) = (Rho(neighbor_p4) - rho0);          // S neighbor
			b_rho2(1) = (Rho(neighbor_p1) - rho0);          // W neighbor
			LU_R24[c].solve(b_rho2, rho_coeff_24);

			b_u_2(0) = (U(neighbor_p4) - u0);   // S neighbor
			b_u_2(1) = (U(neighbor_p1) - u0);   // W neighbor
			LU_R24[c].solve(b_u_2, u_coeff_24);

			b_v_2(0) = (V(neighbor_p4) - v0);   // S neighbor
			b_v_2(1) = (V(neighbor_p1) - v0);   // W neighbor
			LU_R24[c].solve(b_v_2, v_coeff_24);

			b_p_2(0) = (P(neighbor_p4) - p0);               // S neighbor
			b_p_2(1) = (P(neighbor_p1) - p0);               // W neighbor
			LU_R24[c].solve(b_p_2, p_coeff_24);

			b_phi_2(0) = (Phi(neighbor_p4) - phi0);               // S neighbor
			b_phi_2(1) = (Phi(neighbor_p1) - phi0);               // W neighbor
			LU_R24[c].solve(b_phi_2, phi_coeff_24);
 //  			//std::cout<<"interior end"<<std::endl;    
        } // End of interior cell loop 
        
        
        else {

            bool W_face = false, E_face = false, N_face = false, S_face = false;
            					//std::cout<<"c bou: "<<c<<std::endl;
            if (!(is_corner_cell[c])) {
            					//std::cout<<"c non corner: "<<c<<std::endl;
//	   			//std::cout<<"boundary "<<std::endl;    
                
				if(cell->face(0)->at_boundary()) { W_face = true;}
				if(cell->face(1)->at_boundary()) { E_face = true;}
				if(cell->face(2)->at_boundary()) { S_face = true;}
				if(cell->face(3)->at_boundary()) { N_face = true;}
                
                // =====================================================================
                // r = 4 stencil (boundary)
                // =====================================================================
                        	//if(c == 127) //std::cout<<"c non corner 4: 0"<<c<<std::endl;            
                Vector<double> d_rho; Vector<double> d_U; Vector<double> d_V; Vector<double> d_p;
				Vector<double> d_phi;
                
                d_rho.reinit(7); d_U.reinit(7); d_V.reinit(7); d_p.reinit(7);
				d_phi.reinit(7);
                
                index = 0;
                
                if (!W_face) {
                    
                        d_rho(index)   = Rho(neighbor_p1) - rho0;
                        d_U(index) = U(neighbor_p1) - u0;
                        d_V(index) = V(neighbor_p1) - v0;
                        d_p(index)     = P(neighbor_p1) - p0;
                        d_phi(index)     = Phi(neighbor_p1) - phi0;
                        
                        index++; 
                }
                
                if (!N_face) {
                    
                    d_rho(index)   = Rho(neighbor_p2) - rho0;
                    d_U(index) = U(neighbor_p2) - u0;
                    d_V(index) = V(neighbor_p2) - v0;
                    d_p(index)     = P(neighbor_p2) - p0;
                    d_phi(index)     = Phi(neighbor_p2) - phi0;
                    
                    index++; 
                    
                }
                
                if (!E_face) {
                    
                    d_rho(index)   = Rho(neighbor_p3) - rho0;
                    d_U(index) = U(neighbor_p3) - u0;
                    d_V(index) = V(neighbor_p3) - v0;
                    d_p(index)     = P(neighbor_p3) - p0;
                    d_phi(index)     = Phi(neighbor_p3) - phi0;
                    
                    index++; 
                    
                }
                
                if (!S_face) {
                    
                    d_rho(index)   = Rho(neighbor_p4) - rho0;
                    d_U(index) = U(neighbor_p4) - u0;
                    d_V(index) = V(neighbor_p4) - v0;
                    d_p(index)     = P(neighbor_p4) - p0;
                    d_phi(index)     = Phi(neighbor_p4) - phi0;
                    
                    index++; 
                    
                }
                
                d_rho(3)   = 0.0; d_rho(4)   = 0.0; d_rho(5)   = 0.0; d_rho(6)   = 0.0;
                d_U(3) = 0.0; d_U(4) = 0.0; d_U(5) = 0.0; d_U(6) = 0.0;
                d_V(3) = 0.0; d_V(4) = 0.0; d_V(5) = 0.0; d_V(6) = 0.0;
                d_p(3)     = 0.0; d_p(4)     = 0.0; d_p(5)     = 0.0; d_p(6)     = 0.0;
                d_phi(3)     = 0.0; d_phi(4)     = 0.0; d_phi(5)     = 0.0; d_phi(6)     = 0.0;
                            	//if(c == 127) //std::cout<<"c non corner 4: 1"<<c<<std::endl;            
                index = 0;
                
				// Least Squares Part  
			
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
				
				b_rho_4.reinit(ROWS); b_u_4.reinit(ROWS); b_v_4.reinit(ROWS); b_p_4.reinit(ROWS);
				b_phi_4.reinit(ROWS);
				
				// vertex neighbors
            	//if(c == 127) //std::cout<<"c non corner 4: 2"<<c<<std::endl;            				
				for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			

					local_neighbor_dof_indices[0] = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];	

					b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
            	//if(c == 127) //std::cout<<"c non corner 4: 3"<<c<<std::endl;            												
				// neighbors of neighbors 
				
				if (WW) {

					local_neighbor_dof_indices[0] = WW_index;

					b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
            	//if(c == 127) //std::cout<<"c non corner 4: 4"<<c<<std::endl;            								
				if (NN) {

					local_neighbor_dof_indices[0] = NN_index;

					b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
            	//if(c == 127) //std::cout<<"c non corner 4: 5"<<c<<std::endl;            								
				if (EE) {
            	//if(c == 127) //std::cout<<"c non corner 4: 5: "<<EE_index<<std::endl;            								
					local_neighbor_dof_indices[0] = EE_index;
            	//if(c == 127) //std::cout<<"c non corner 4: 5: "<<EE_index<<std::endl;            								
					b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
            	//if(c == 127) //std::cout<<"c non corner 4: 6"<<c<<std::endl;            								
				if (SS) {

					local_neighbor_dof_indices[0] = SS_index;

					b_rho_4(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
					b_u_4(index) = U(local_neighbor_dof_indices[0] ) - u0;
					b_v_4(index) = V(local_neighbor_dof_indices[0] ) - v0;
					b_p_4(index)     = P(local_neighbor_dof_indices[0] ) - p0;
					b_phi_4(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
					index++; 
				}
                            	//if(c == 127) //std::cout<<"c non corner 4: 7"<<c<<std::endl;            
                
                CLS_R4[c].solve(b_rho_4, d_rho, rho_coeff_4);
                CLS_R4[c].solve(b_u_4, d_U, u_coeff_4);
                CLS_R4[c].solve(b_v_4, d_V, v_coeff_4);
                CLS_R4[c].solve(b_p_4, d_p, p_coeff_4);
                CLS_R4[c].solve(b_phi_4, d_phi, phi_coeff_4);
            	//if(c == 127) //std::cout<<"c non corner 4: "<<c<<std::endl;
				// =====================================================================
                // r = 3 center stencil (boundary)
                // =====================================================================
                
                b_rho_3.reinit(5); b_u_3.reinit(5); b_v_3.reinit(5); b_p_3.reinit(5);  
				b_phi_3.reinit(5);  
                
                index = 0; 
                
                if (!W_face) {
                    
                    // P1 neighbor
                    
                    b_rho_3(index)    = (Rho(neighbor_p1) - rho0);                 
                    b_u_3(index)  = (U(neighbor_p1) - u0);                
                    b_v_3(index)  = (V(neighbor_p1) - v0);                
                    b_p_3(index)      = (P(neighbor_p1) - p0);                
                    b_phi_3(index)      = (Phi(neighbor_p1) - phi0);
                    
                    index++; 
                }
                
                if (!N_face) {

                    // P2 neighbor
                    
                    b_rho_3(index)    = (Rho(neighbor_p2) - rho0);                 
                    b_u_3(index)  = (U(neighbor_p2) - u0);                
                    b_v_3(index)  = (V(neighbor_p2) - v0);                
                    b_p_3(index)      = (P(neighbor_p2) - p0);                
                    b_phi_3(index)      = (Phi(neighbor_p2) - phi0);                
                    
                    index++; 
                
                }
                
                if (!E_face) {

                    // P3 neighbor
                    
                    b_rho_3(index)    = (Rho(neighbor_p3) - rho0);                 
                    b_u_3(index)  = (U(neighbor_p3) - u0);                
                    b_v_3(index)  = (V(neighbor_p3) - v0);                
                    b_p_3(index)      = (P(neighbor_p3) - p0);                
                    b_phi_3(index)      = (Phi(neighbor_p3) - phi0);                
                    
                    index++; 
                }
                
                if (!S_face) {

                    // P4 neighbor
                    
                    b_rho_3(index)    = (Rho(neighbor_p4) - rho0);                 
                    b_u_3(index)  = (U(neighbor_p4) - u0);                
                    b_v_3(index)  = (V(neighbor_p4) - v0);                
                    b_p_3(index)      = (P(neighbor_p4) - p0);                
                    b_phi_3(index)      = (Phi(neighbor_p4) - phi0);                
                    
                    index++; 
                }
                
                
                // Transmissive boundary conditions 
                
				b_rho_3(3) = 0.0; b_u_3(3) = 0.0; b_v_3(3) = 0.0; b_p_3(3) = 0.0;
				b_phi_3(3) = 0.0;
				b_rho_3(4) = 0.0; b_u_3(4) = 0.0; b_v_3(4) = 0.0; b_p_3(4) = 0.0;
				b_phi_3(4) = 0.0;
                
            
				CLS_R3[c].solve(b_rho_3, rho_coeff_3); 
				CLS_R3[c].solve(b_u_3, u_coeff_3);
				CLS_R3[c].solve(b_v_3, v_coeff_3);
				CLS_R3[c].solve(b_p_3, p_coeff_3);
				CLS_R3[c].solve(b_phi_3, phi_coeff_3);
            	//if(c == 127) //std::cout<<"c non corner 3: "<<c<<std::endl;
				// =====================================================================
                // r = 3 stencil 1 (boundary)
                // =====================================================================
                
                if (is_admissible_R31[c]) {
                    
                    index = 0; 
                    if ( N_face || S_face ) {                    
                        d_rho_31.reinit(5); d_u_31.reinit(5); d_v_31.reinit(5); d_p_31.reinit(5);  
						d_phi_31.reinit(5);  
                    }
            
                    else {
                        d_rho_31.reinit(4); d_u_31.reinit(4); d_v_31.reinit(4); d_p_31.reinit(4);
						d_phi_31.reinit(4);
                    }
                    
                    if (!W_face) {
                    
                        d_rho_31(index)  = Rho(neighbor_p1) - rho0;                // P1 neighbor 
                        d_u_31(index)  = U(neighbor_p1) - u0;          // P1 neighbor 
                        d_v_31(index)  = V(neighbor_p1) - v0;          // P1 neighbor 
                        d_p_31(index)  = P(neighbor_p1) - p0;                      // P1 neighbor 
                        d_phi_31(index)  = Phi(neighbor_p1) - phi0;                      // P1 neighbor 
                        
                        index++;
                    }
                    
                    if (!N_face) {
                    
                        d_rho_31(index)  = Rho(neighbor_p2) - rho0;                // P2 neighbor 
                        d_u_31(index)  = U(neighbor_p2) - u0;          // P2 neighbor 
                        d_v_31(index)  = V(neighbor_p2) - v0;          // P2 neighbor 
                        d_p_31(index)  = P(neighbor_p2) - p0;                      // P2 neighbor 
                        d_phi_31(index)  = Phi(neighbor_p2) - phi0;                      // P2 neighbor 
                        
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_rho_31(index)  = Rho(neighbor_p4) - rho0;                // P4 neighbor 
                        d_u_31(index)  = U(neighbor_p4) - u0;          // P4 neighbor 
                        d_v_31(index)  = V(neighbor_p4) - v0;          // P4 neighbor 
                        d_p_31(index)  = P(neighbor_p4) - p0;                      // P4 neighbor 
                        d_phi_31(index)  = Phi(neighbor_p4) - phi0;                      // P4 neighbor 
                        
                        index++;
                    }
                    
                    if (WW) {

						local_neighbor_dof_indices[0] = WW_index;
                    
                        d_rho_31(index)    = Rho(local_neighbor_dof_indices[0] ) - rho0;              // S1 neighbor 
                        d_u_31(index)  = U(local_neighbor_dof_indices[0] ) - u0;          // S1 neighbor 
                        d_v_31(index)  = V(local_neighbor_dof_indices[0] ) - v0;          // S1 neighbor 
                        d_p_31(index)      = P(local_neighbor_dof_indices[0] ) - p0;                  // S1 neighbor 
                        d_phi_31(index)      = Phi(local_neighbor_dof_indices[0] ) - phi0;                  // S1 neighbor 
                        
                        index++;
                    }
                    if ( S_face || N_face ) {                    
                        
                        // Transmissive boundary conditions 
                        
                        d_rho_31(index)    = 0.0;   d_rho_31(index+1)    = 0.0;            
                        d_u_31(index)  = 0.0;   d_u_31(index+1)  = 0.0;      
                        d_v_31(index)  = 0.0;   d_v_31(index+1)  = 0.0;    
                        d_p_31(index)      = 0.0;   d_p_31(index+1)  = 0.0;
                        d_phi_31(index)      = 0.0;   d_phi_31(index+1)  = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][0].size(); 
				
					b_rho_31.reinit(ROWS); b_u_31.reinit(ROWS); b_v_31.reinit(ROWS); b_p_31.reinit(ROWS);
					b_phi_31.reinit(ROWS);

					index = 0; 
					
					// vertex neighbor of cell at face 0 
					for (unsigned int d = 0; d < cell_neighbor_index[c][0].size(); ++d) {						

						local_neighbor_dof_indices[0] = global_to_local_index_map[cell_neighbor_index[c][0][d] ];

						b_rho_31(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
						b_u_31(index) = U(local_neighbor_dof_indices[0] ) - u0;
						b_v_31(index) = V(local_neighbor_dof_indices[0] ) - v0;
						b_p_31(index)     = P(local_neighbor_dof_indices[0] ) - p0;
						b_phi_31(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
						index++; 
					}

					CLS_R31[c].solve(b_rho_31, d_rho_31, rho_coeff_31);
					CLS_R31[c].solve(b_u_31, d_u_31, u_coeff_31);
					CLS_R31[c].solve(b_v_31, d_v_31, v_coeff_31);
					CLS_R31[c].solve(b_p_31, d_p_31, p_coeff_31);
					CLS_R31[c].solve(b_phi_31, d_phi_31, phi_coeff_31);
				
				}
                
                else {
					rho_coeff_31 = rho_coeff_3; 
					u_coeff_31 = u_coeff_3;
					v_coeff_31 = v_coeff_3;
					p_coeff_31 = p_coeff_3;
					phi_coeff_31 = phi_coeff_3;
                }
            	//if(c == 127) //std::cout<<"c non corner 31: "<<c<<std::endl;                
				// =====================================================================
                // r = 3 stencil 2 (boundary)
                // =====================================================================
                
                if (is_admissible_R32[c]) {
                    
                    index = 0; 
                    if ( W_face || E_face ) {                    
                        d_rho_32.reinit(5); d_u_32.reinit(5); d_v_32.reinit(5); d_p_32.reinit(5);  
						d_phi_32.reinit(5);  
                    }
            
                    else {
                        d_rho_32.reinit(4); d_u_32.reinit(4); d_v_32.reinit(4); d_p_32.reinit(4);
						d_phi_32.reinit(4);
                    }
                    
                    if (!W_face) {
                    
                        d_rho_32(index)  = Rho(neighbor_p1) - rho0;                // W neighbor 
                        d_u_32(index)  = U(neighbor_p1) - u0;          // W neighbor 
                        d_v_32(index)  = V(neighbor_p1) - v0;          // W neighbor 
                        d_p_32(index)  = P(neighbor_p1) - p0;                      // W neighbor 
                        d_phi_32(index)  = Phi(neighbor_p1) - phi0;                      // W neighbor 
                        
                        index++;
                    }
                    
                    if (!N_face) {
                    
                        d_rho_32(index)  = Rho(neighbor_p2) - rho0;                // N neighbor 
                        d_u_32(index)  = U(neighbor_p2) - u0;          // N neighbor 
                        d_v_32(index)  = V(neighbor_p2) - v0;          // N neighbor 
                        d_p_32(index)  = P(neighbor_p2) - p0;                      // N neighbor 
                        d_phi_32(index)  = Phi(neighbor_p2) - phi0;                      // N neighbor 
                        
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_rho_32(index)  = Rho(neighbor_p3) - rho0;                // P3 neighbor 
                        d_u_32(index)  = U(neighbor_p3) - u0;          // P3 neighbor 
                        d_v_32(index)  = V(neighbor_p3) - v0;          // P3 neighbor 
                        d_p_32(index)  = P(neighbor_p3) - p0;                      // P3 neighbor 
                        d_phi_32(index)  = Phi(neighbor_p3) - phi0;                      // P3 neighbor 
                        
                        index++;
                    }
                    
                    if (NN) {

						local_neighbor_dof_indices[0] = NN_index;
                     
                        d_rho_32(index)    = Rho(local_neighbor_dof_indices[0] ) - rho0;              // S3 neighbor 
                        d_u_32(index)  = U(local_neighbor_dof_indices[0] ) - u0;          // S3 neighbor 
                        d_v_32(index)  = V(local_neighbor_dof_indices[0] ) - v0;          // S3 neighbor 
                        d_p_32(index)      = P(local_neighbor_dof_indices[0] ) - p0;                  // S3 neighbor 
                        d_phi_32(index)      = Phi(local_neighbor_dof_indices[0] ) - phi0;                  // S3 neighbor 
                        
                        index++;
                    }
                    
                    if ( W_face || E_face ) {
                        
                        // Transmissive boundary conditions 
                        
                        d_rho_32(index)    = 0.0;   d_rho_32(index+1)    = 0.0;            
                        d_u_32(index)  = 0.0;   d_u_32(index+1)  = 0.0;      
                        d_v_32(index)  = 0.0;   d_v_32(index+1)  = 0.0;    
                        d_p_32(index)      = 0.0;   d_p_32(index+1)  = 0.0;
                        d_phi_32(index)      = 0.0;   d_phi_32(index+1)  = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][3].size(); 

					b_rho_32.reinit(ROWS); b_u_32.reinit(ROWS); b_v_32.reinit(ROWS); b_p_32.reinit(ROWS);
					b_phi_32.reinit(ROWS);

					index = 0; 		

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][3].size(); ++d) {						

						local_neighbor_dof_indices[0] = global_to_local_index_map[cell_neighbor_index[c][3][d] ];

						b_rho_32(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
						b_u_32(index) = U(local_neighbor_dof_indices[0] ) - u0;
						b_v_32(index) = V(local_neighbor_dof_indices[0] ) - v0;
						b_p_32(index)     = P(local_neighbor_dof_indices[0] ) - p0;
						b_phi_32(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
						index++; 
					}
                    
					CLS_R32[c].solve(b_rho_32, d_rho_32, rho_coeff_32);
					CLS_R32[c].solve(b_u_32, d_u_32, u_coeff_32);
					CLS_R32[c].solve(b_v_32, d_v_32, v_coeff_32);
					CLS_R32[c].solve(b_p_32, d_p_32, p_coeff_32); 
					CLS_R32[c].solve(b_phi_32, d_phi_32, phi_coeff_32); 
					
				}
                
                else {
					rho_coeff_32 = rho_coeff_3; 
					u_coeff_32 = u_coeff_3;
					v_coeff_32 = v_coeff_3;
					p_coeff_32 = p_coeff_3;
					phi_coeff_32 = phi_coeff_3;
				}
            	//if(c == 127) //std::cout<<"c non corner 32: "<<c<<std::endl;                                
                // =====================================================================
                // r = 3 stencil 3 (boundary)
                // =====================================================================
                
                if (is_admissible_R33[c]) {
                    
                    index = 0; 
                    if ( W_face || E_face ) {                    
                        d_rho_33.reinit(5); d_u_33.reinit(5); d_v_33.reinit(5); d_p_33.reinit(5);  
						d_phi_33.reinit(5);  
                    }
            
                    else {
                        d_rho_33.reinit(4); d_u_33.reinit(4); d_v_33.reinit(4); d_p_33.reinit(4);
						d_phi_33.reinit(4);
                    }
                    
                    if (!W_face) {
                    
                        d_rho_33(index)  = Rho(neighbor_p1) - rho0;                // W neighbor 
                        d_u_33(index)  = U(neighbor_p1) - u0;          // W neighbor 
                        d_v_33(index)  = V(neighbor_p1) - v0;          // W neighbor 
                        d_p_33(index)  = P(neighbor_p1) - p0;                      // W neighbor 
                        d_phi_33(index)  = Phi(neighbor_p1) - phi0;                      // W neighbor 
                        
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_rho_33(index)  = Rho(neighbor_p3) - rho0;                // E neighbor 
                        d_u_33(index)  = U(neighbor_p3) - u0;          // E neighbor 
                        d_v_33(index)  = V(neighbor_p3) - v0;          // E neighbor 
                        d_p_33(index)  = P(neighbor_p3) - p0;                      // E neighbor 
                        d_phi_33(index)  = Phi(neighbor_p3) - phi0;                      // E neighbor 
                        
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_rho_33(index)  = Rho(neighbor_p4) - rho0;                // S neighbor 
                        d_u_33(index)  = U(neighbor_p4) - u0;          // S neighbor 
                        d_v_33(index)  = V(neighbor_p4) - v0;          // S neighbor 
                        d_p_33(index)  = P(neighbor_p4) - p0;                      // S neighbor 
                        d_phi_33(index)  = Phi(neighbor_p4) - phi0;                      // S neighbor 
                        
                        index++;
                    }
                    
                    if (SS) {

						local_neighbor_dof_indices[0] = SS_index;
                    
                        d_rho_33(index)    = Rho(local_neighbor_dof_indices[0] ) - rho0;              // SS neighbor 
                        d_u_33(index)  = U(local_neighbor_dof_indices[0] ) - u0;          // SS neighbor 
                        d_v_33(index)  = V(local_neighbor_dof_indices[0] ) - v0;          // SS neighbor 
                        d_p_33(index)      = P(local_neighbor_dof_indices[0] ) - p0;                  // SS neighbor 
                        d_phi_33(index)      = Phi(local_neighbor_dof_indices[0] ) - phi0;                  // SS neighbor 
                        
                        index++;
                    }
                    
                    if ( W_face || E_face ) {                    
                        
                        // Transmissive boundary conditions 
                        
                        d_rho_33(index)    = 0.0;   d_rho_33(index+1)    = 0.0;            
                        d_u_33(index)  = 0.0;   d_u_33(index+1)  = 0.0;      
                        d_v_33(index)  = 0.0;   d_v_33(index+1)  = 0.0;    
                        d_p_33(index)      = 0.0;   d_p_33(index+1)  = 0.0;
                        d_phi_33(index)      = 0.0;   d_phi_33(index+1)  = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][2].size(); 

					b_rho_33.reinit(ROWS); b_u_33.reinit(ROWS); b_v_33.reinit(ROWS); b_p_33.reinit(ROWS);
					b_phi_33.reinit(ROWS);

					index = 0; 							

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][2].size(); ++d) {						

						local_neighbor_dof_indices[0] = global_to_local_index_map[cell_neighbor_index[c][2][d] ];

						b_rho_33(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
						b_u_33(index) = U(local_neighbor_dof_indices[0] ) - u0;
						b_v_33(index) = V(local_neighbor_dof_indices[0] ) - v0;
						b_p_33(index)     = P(local_neighbor_dof_indices[0] ) - p0;
						b_phi_33(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
						index++; 
					}

					CLS_R33[c].solve(b_rho_33, d_rho_33, rho_coeff_33);
					CLS_R33[c].solve(b_u_33, d_u_33, u_coeff_33);
					CLS_R33[c].solve(b_v_33, d_v_33, v_coeff_33);
					CLS_R33[c].solve(b_p_33, d_p_33, p_coeff_33); 
					CLS_R33[c].solve(b_phi_33, d_phi_33, phi_coeff_33); 
				}
                
                else {
					rho_coeff_33 = rho_coeff_3; 
					u_coeff_33 = u_coeff_3;
					v_coeff_33 = v_coeff_3;
					p_coeff_33 = p_coeff_3;
					phi_coeff_33 = phi_coeff_3;
                }
                            	//if(c == 127) //std::cout<<"c non corner 33: "<<c<<std::endl;                
                // =====================================================================
                // r = 3 stencil 4 (boundary)
                // =====================================================================
                
                if (is_admissible_R34[c]) {
                    
                    index = 0; 
                    if ( N_face || S_face ) {                                        
                        d_rho_34.reinit(5); d_u_34.reinit(5); d_v_34.reinit(5); d_p_34.reinit(5);  
						d_phi_34.reinit(5);  
                    }
            
                    else {
                        d_rho_34.reinit(4); d_u_34.reinit(4); d_v_34.reinit(4); d_p_34.reinit(4);
						d_phi_34.reinit(4);
                    }
                    
                    if (!N_face) {
                    
                        d_rho_34(index)    = Rho(neighbor_p2) - rho0;              // N neighbor 
                        d_u_34(index)  = U(neighbor_p2) - u0;          // N neighbor 
                        d_v_34(index)  = V(neighbor_p2) - v0;          // N neighbor 
                        d_p_34(index)      = P(neighbor_p2) - p0;                  // N neighbor 
                        d_phi_34(index)      = Phi(neighbor_p2) - phi0;                  // N neighbor 
                        
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_rho_34(index)    = Rho(neighbor_p3) - rho0;              // E neighbor 
                        d_u_34(index)  = U(neighbor_p3) - u0;          // E neighbor 
                        d_v_34(index)  = V(neighbor_p3) - v0;          // E neighbor 
                        d_p_34(index)      = P(neighbor_p3) - p0;                  // E neighbor 
                        d_phi_34(index)      = Phi(neighbor_p3) - phi0;                  // E neighbor 
                        
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_rho_34(index)    = Rho(neighbor_p4) - rho0;              // S neighbor 
                        d_u_34(index)  = U(neighbor_p4) - u0;          // S neighbor 
                        d_v_34(index)  = V(neighbor_p4) - v0;          // S neighbor 
                        d_p_34(index)      = P(neighbor_p4) - p0;                  // S neighbor 
                        d_phi_34(index)      = Phi(neighbor_p4) - phi0;                  // S neighbor 
                        
                        index++;
                    }
                    
                    if (EE) {

						local_neighbor_dof_indices[0] = EE_index;
                    
                        d_rho_34(index)    = Rho(local_neighbor_dof_indices[0] ) - rho0;              // S6 neighbor 
                        d_u_34(index)  = U(local_neighbor_dof_indices[0] ) - u0;          // S6 neighbor 
                        d_v_34(index)  = V(local_neighbor_dof_indices[0] ) - v0;          // S6 neighbor 
                        d_p_34(index)      = P(local_neighbor_dof_indices[0] ) - p0;                  // S6 neighbor 
                        d_phi_34(index)      = Phi(local_neighbor_dof_indices[0] ) - phi0;                  // S6 neighbor 
                        
                        index++;
                    }
                    
                    if ( N_face || S_face ) {                                        
                        
                        // Transmissive boundary conditions 
                        
                        d_rho_34(index)    = 0.0;   d_rho_34(index+1)    = 0.0;            
                        d_u_34(index)  = 0.0;   d_u_34(index+1)  = 0.0;      
                        d_v_34(index)  = 0.0;   d_v_34(index+1)  = 0.0;    
                        d_p_34(index)      = 0.0;   d_p_34(index+1)  = 0.0;
                        d_phi_34(index)      = 0.0;   d_phi_34(index+1)  = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][1].size(); 

					b_rho_34.reinit(ROWS); b_u_34.reinit(ROWS); b_v_34.reinit(ROWS); b_p_34.reinit(ROWS);
					b_phi_34.reinit(ROWS);

					index = 0; 		

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][1].size(); ++d) {						

						local_neighbor_dof_indices[0] = global_to_local_index_map[cell_neighbor_index[c][1][d] ];

						b_rho_34(index)   = Rho(local_neighbor_dof_indices[0] ) - rho0; 
						b_u_34(index) = U(local_neighbor_dof_indices[0] ) - u0;
						b_v_34(index) = V(local_neighbor_dof_indices[0] ) - v0;
						b_p_34(index)     = P(local_neighbor_dof_indices[0] ) - p0;
						b_phi_34(index)     = Phi(local_neighbor_dof_indices[0] ) - phi0;
						index++; 
					}                   
					
                    CLS_R34[c].solve(b_rho_34, d_rho_34, rho_coeff_34);
                    CLS_R34[c].solve(b_u_34, d_u_34, u_coeff_34);
                    CLS_R34[c].solve(b_v_34, d_v_34, v_coeff_34);
                    CLS_R34[c].solve(b_p_34, d_p_34, p_coeff_34); 
                    CLS_R34[c].solve(b_phi_34, d_phi_34, phi_coeff_34); 
                }
                
                else {
					
					rho_coeff_34 = rho_coeff_3; 
					u_coeff_34 = u_coeff_3;
					v_coeff_34 = v_coeff_3;
					p_coeff_34 = p_coeff_3;
					phi_coeff_34 = phi_coeff_3;
                }
            	//if(c == 127) //std::cout<<"c non corner 34: "<<c<<std::endl;                              
                // =====================================================================
                // r = 2 stencil 1 (boundary)
                // =====================================================================
                
                if (W_face) { // WEST boundary 
                
                    b_rho2(0) = (Rho(neighbor_p2) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_u_2(0) = (U(neighbor_p2) - u0);   // P2 neighbor
                    b_u_2(1) = 0.0; 
                    
                    b_v_2(0) = (V(neighbor_p2) - v0);   // P2 neighbor
                    b_v_2(1) = 0.0; 
                    
                    b_p_2(0) = (P(neighbor_p2) - p0);   // P2 neighbor
                    b_p_2(1) = 0.0; 

                    b_phi_2(0) = (Phi(neighbor_p2) - phi0);   // P2 neighbor
                    b_phi_2(1) = 0.0; 
                }
                
                else if (N_face) { // NORTH boundary 
                
                    b_rho2(0) = (Rho(neighbor_p1) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_u_2(0) = (U(neighbor_p1) - u0);   // P2 neighbor
                    b_u_2(1) = 0.0; 
                    
                    b_v_2(0) = (V(neighbor_p1) - v0);   // P2 neighbor
                    b_v_2(1) = 0.0; 
                    
                    b_p_2(0) = (P(neighbor_p1) - p0);   // P2 neighbor
                    b_p_2(1) = 0.0; 

                    b_phi_2(0) = (Phi(neighbor_p1) - phi0);   // P2 neighbor
                    b_phi_2(1) = 0.0; 
                }
                
                else {
                    
                    b_rho2(0)  = (Rho(neighbor_p1) - rho0);   // P1 neighbor 
                    b_rho2(1)  = (Rho(neighbor_p2) - rho0);   // P2 neighbor

                    b_u_2(0)  = (U(neighbor_p1) - u0);   // P1 neighbor 
                    b_u_2(1)  = (U(neighbor_p2) - u0);   // P2 neighbor

                    b_v_2(0)  = (V(neighbor_p1) - v0);   // P1 neighbor 
                    b_v_2(1)  = (V(neighbor_p2) - v0);   // P2 neighbor

                    b_p_2(0)  = (P(neighbor_p1) - p0);   // P1 neighbor 
                    b_p_2(1)  = (P(neighbor_p2) - p0);   // P2 neighbor

                    b_phi_2(0)  = (Phi(neighbor_p1) - phi0);   // P1 neighbor 
                    b_phi_2(1)  = (Phi(neighbor_p2) - phi0);   // P2 neighbor

                }
                
				LU_R21[c].solve(b_rho2, rho_coeff_21);
                LU_R21[c].solve(b_u_2, u_coeff_21);
                LU_R21[c].solve(b_v_2, v_coeff_21);
                LU_R21[c].solve(b_p_2, p_coeff_21);
                LU_R21[c].solve(b_phi_2, phi_coeff_21);
            	//if(c == 127) //std::cout<<"c non corner 21: "<<c<<std::endl;                
                // =====================================================================
                // r = 2 stencil 2 (boundary)
                // =====================================================================
                
                if (E_face) { // EAST boundary 
                
                    b_rho2(0) = (Rho(neighbor_p2) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_u_2(0) = (U(neighbor_p2) - u0);   // P2 neighbor
                    b_u_2(1) = 0.0; 
                    
                    b_v_2(0) = (V(neighbor_p2) - v0);   // P2 neighbor
                    b_v_2(1) = 0.0; 
                    
                    b_p_2(0) = (P(neighbor_p2) - p0);   // P2 neighbor
                    b_p_2(1) = 0.0; 

                    b_phi_2(0) = (Phi(neighbor_p2) - phi0);   // P2 neighbor
                    b_phi_2(1) = 0.0; 
                }
                
                else if (N_face) { // NORTH boundary 
                
                    b_rho2(0) = (Rho(neighbor_p3) - rho0);   // P3 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_u_2(0) = (U(neighbor_p3) - u0);   // P3 neighbor
                    b_u_2(1) = 0.0; 
                    
                    b_v_2(0) = (V(neighbor_p3) - v0);   // P3 neighbor
                    b_v_2(1) = 0.0; 
                    
                    b_p_2(0) = (P(neighbor_p3) - p0);   // P2 neighbor
                    b_p_2(1) = 0.0; 

                    b_phi_2(0) = (Phi(neighbor_p3) - phi0);   // P2 neighbor
                    b_phi_2(1) = 0.0; 

                }
                
                else {
                    
                    b_rho2(0) = (Rho(neighbor_p2) - rho0);   // P2 neighbor
                    b_rho2(1) = (Rho(neighbor_p3) - rho0);   // P3 neighbor
                    
                    b_u_2(0) = (U(neighbor_p2) - u0);   // P2 neighbor
                    b_u_2(1) = (U(neighbor_p3) - u0);   // P3 neighbor

                    b_v_2(0) = (V(neighbor_p2) - v0);   // P2 neighbor
                    b_v_2(1) = (V(neighbor_p3) - v0);   // P3 neighbor

                    b_p_2(0) = (P(neighbor_p2) - p0);   // P2 neighbor
                    b_p_2(1) = (P(neighbor_p3) - p0);   // P3 neighbor

                    b_phi_2(0) = (Phi(neighbor_p2) - phi0);   // P2 neighbor
                    b_phi_2(1) = (Phi(neighbor_p3) - phi0);   // P3 neighbor

                }
                
                LU_R22[c].solve(b_rho2, rho_coeff_22);
                LU_R22[c].solve(b_u_2, u_coeff_22);
                LU_R22[c].solve(b_v_2, v_coeff_22);
                LU_R22[c].solve(b_p_2, p_coeff_22);
                LU_R22[c].solve(b_phi_2, phi_coeff_22);
            	//if(c == 127) //std::cout<<"c non corner 22: "<<c<<std::endl;                                

                // =====================================================================
                // r = 2 stencil 3 (boundary)
                // =====================================================================
                
                if (E_face) { // EAST boundary 
                
                    b_rho2(0) = (Rho(neighbor_p4) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_u_2(0) = (U(neighbor_p4) - u0);   // P2 neighbor
                    b_u_2(1) = 0.0; 
                    
                    b_v_2(0) = (V(neighbor_p4) - v0);   // P2 neighbor
                    b_v_2(1) = 0.0; 
                    
                    b_p_2(0) = (P(neighbor_p4) - p0);   // P2 neighbor
                    b_p_2(1) = 0.0; 

                    b_phi_2(0) = (Phi(neighbor_p4) - phi0);   // P2 neighbor
                    b_phi_2(1) = 0.0; 

                }
                
                else if (S_face) { // SOUTH boundary 
                
                    b_rho2(0) = (Rho(neighbor_p3) - rho0);   // P3 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_u_2(0) = (U(neighbor_p3) - u0);   // P3 neighbor
                    b_u_2(1) = 0.0; 
                    
                    b_v_2(0) = (V(neighbor_p3) - v0);   // P3 neighbor
                    b_v_2(1) = 0.0; 
                    
                    b_p_2(0) = (P(neighbor_p3) - p0);   // P2 neighbor
                    b_p_2(1) = 0.0; 

                    b_phi_2(0) = (Phi(neighbor_p3) - phi0);   // P2 neighbor
                    b_phi_2(1) = 0.0; 

                }
                
                else {
                
                    b_rho2(0) = (Rho(neighbor_p3) - rho0);   // P3 neighbor
                    b_rho2(1) = (Rho(neighbor_p4) - rho0);   // P4 neighbor                                            // P4 neighbor
                     
                    b_u_2(0) = (U(neighbor_p3) - u0);   // P3 neighbor
                    b_u_2(1) = (U(neighbor_p4) - u0);   // P4 neighbor                                            // P4 neighbor
                    
                    b_v_2(0) = (V(neighbor_p3) - v0);   // P3 neighbor
                    b_v_2(1) = (V(neighbor_p4) - v0);   // P4 neighbor                                            // P4 neighbor
                    
                    b_p_2(0) = (P(neighbor_p3) - p0);   // P3 neighbor
                    b_p_2(1) = (P(neighbor_p4) - p0);   // P4 neighbor 

                    b_phi_2(0) = (Phi(neighbor_p3) - phi0);   // P3 neighbor
                    b_phi_2(1) = (Phi(neighbor_p4) - phi0);   // P4 neighbor                                            // P4 neighbor

                }
                
                LU_R23[c].solve(b_rho2, rho_coeff_23);
                LU_R23[c].solve(b_u_2, u_coeff_23);
                LU_R23[c].solve(b_v_2, v_coeff_23);
                LU_R23[c].solve(b_p_2, p_coeff_23);
                LU_R23[c].solve(b_phi_2, phi_coeff_23);
            	//if(c == 127) //std::cout<<"c non corner 23: "<<c<<std::endl;                                                
                // =====================================================================
                // r = 2 stencil 4 (boundary)
                // =====================================================================
                
                if (W_face) { // WEST boundary 
                
                    b_rho2(0) = (Rho(neighbor_p4) - rho0);   // P4 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_u_2(0) = (U(neighbor_p4) - u0);   // P4 neighbor
                    b_u_2(1) = 0.0; 
                    
                    b_v_2(0) = (V(neighbor_p4) - v0);   // P4 neighbor
                    b_v_2(1) = 0.0; 
                    
                    b_p_2(0) = (P(neighbor_p4) - p0);   // P4 neighbor
                    b_p_2(1) = 0.0; 

                    b_phi_2(0) = (Phi(neighbor_p4) - phi0);   // P4 neighbor
                    b_phi_2(1) = 0.0; 

                    
                }
                
                else if (S_face) { // SOUTH boundary 
                
                    b_rho2(0) = (Rho(neighbor_p1) - rho0);   // P1 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_u_2(0) = (U(neighbor_p1) - u0);   // P1 neighbor
                    b_u_2(1) = 0.0; 
                    
                    b_v_2(0) = (V(neighbor_p1) - v0);   // P1 neighbor
                    b_v_2(1) = 0.0; 
                    
                    b_p_2(0) = (P(neighbor_p1) - p0);   // P1 1eighbor
                    b_p_2(1) = 0.0; 

                    b_phi_2(0) = (Phi(neighbor_p1) - phi0);   // P1 1eighbor
                    b_phi_2(1) = 0.0; 


                }
                
                else {
                    b_rho2(0) = (Rho(neighbor_p4) - rho0);   // P4 neighbor
                    b_rho2(1) = (Rho(neighbor_p1) - rho0);   // P1 neighbor

                    b_u_2(0) = (U(neighbor_p4) - u0);   // P4 neighbor
                    b_u_2(1) = (U(neighbor_p1) - u0);   // P1 neighbor
                    
                    b_v_2(0) = (V(neighbor_p4) - v0);   // P4 neighbor
                    b_v_2(1) = (V(neighbor_p1) - v0);   // P1 neighbor
                    
                    b_p_2(0) = (P(neighbor_p4) - p0);   // P4 neighbor
                    b_p_2(1) = (P(neighbor_p1) - p0);   // P1 neighbor

                    b_phi_2(0) = (Phi(neighbor_p4) - phi0);   // P4 neighbor
                    b_phi_2(1) = (Phi(neighbor_p1) - phi0);   // P1 neighbor
                    
                }
                
                LU_R24[c].solve(b_rho2, rho_coeff_24);
                LU_R24[c].solve(b_u_2, u_coeff_24);
                LU_R24[c].solve(b_v_2, v_coeff_24);
                LU_R24[c].solve(b_p_2, p_coeff_24);
                LU_R24[c].solve(b_phi_2, phi_coeff_24);
            	//if(c == 127) //std::cout<<"c non corner 24: "<<c<<std::endl;                                
//	   			//std::cout<<"boundary end"<<std::endl;                    
            } // Non-corner boundary cell loop 
            
            else {
                
                // Corner boundary cells - reduce to first order 
                
                rho_coeff_4 = 0.0;
				rho_coeff_3 = 0.0;
				rho_coeff_31 = 0.0;
				rho_coeff_32 = 0.0;
                rho_coeff_33 = 0.0;  
                rho_coeff_34 = 0.0;  
				rho_coeff_21 = 0.0;
				rho_coeff_22 = 0.0;
				rho_coeff_23 = 0.0;
				rho_coeff_24 = 0.0;
				
				u_coeff_4 = 0.0;
				u_coeff_3 = 0.0;
				u_coeff_31 = 0.0;
				u_coeff_32 = 0.0;
                u_coeff_32 = 0.0;  
                u_coeff_34 = 0.0;  
				u_coeff_21 = 0.0;
				u_coeff_22 = 0.0;
				u_coeff_23 = 0.0;
				u_coeff_24 = 0.0;
				
				v_coeff_4 = 0.0;
				v_coeff_3 = 0.0;
				v_coeff_31 = 0.0;
				v_coeff_32 = 0.0;
                v_coeff_33 = 0.0;  
                v_coeff_34 = 0.0;  
				v_coeff_21 = 0.0;
				v_coeff_22 = 0.0;
				v_coeff_23 = 0.0;
				v_coeff_24 = 0.0;
				
				p_coeff_4 = 0.0;
				p_coeff_3 = 0.0;
				p_coeff_31 = 0.0;
				p_coeff_32 = 0.0;
                p_coeff_33 = 0.0;  
                p_coeff_34 = 0.0;  
				p_coeff_21 = 0.0;
				p_coeff_22 = 0.0;
				p_coeff_23 = 0.0;
				p_coeff_24 = 0.0;

				phi_coeff_4 = 0.0;
				phi_coeff_3 = 0.0;
				phi_coeff_31 = 0.0;
				phi_coeff_32 = 0.0;
                phi_coeff_33 = 0.0;  
                phi_coeff_34 = 0.0;  
				phi_coeff_21 = 0.0;
				phi_coeff_22 = 0.0;
				phi_coeff_23 = 0.0;
				phi_coeff_24 = 0.0;
                
            } // Corner boundary cell loop 
            
        } // End of boundary cell loop 
        
        // Find the smoothness indicators 
		// =====================================================================
		// r = 4 Stencil 
		// =====================================================================
        
		IS_RHO(0)   = compute_fourth_order_smoothness_indicator(rho_coeff_4);
		IS_U(0) = compute_fourth_order_smoothness_indicator(u_coeff_4);
		IS_V(0) = compute_fourth_order_smoothness_indicator(v_coeff_4);
		IS_P(0)     = compute_fourth_order_smoothness_indicator(p_coeff_4);
		IS_PHI(0)     = compute_fourth_order_smoothness_indicator(phi_coeff_4);
		// =====================================================================
		// r = 3
		// =====================================================================
        
		IS_RHO(1)   = compute_third_order_smoothness_indicator(rho_coeff_3); 
		IS_U(1) = compute_third_order_smoothness_indicator(u_coeff_3);
		IS_V(1) = compute_third_order_smoothness_indicator(v_coeff_3);
		IS_P(1)     = compute_third_order_smoothness_indicator(p_coeff_3);
		IS_PHI(1)     = compute_third_order_smoothness_indicator(phi_coeff_3);
        
		// =====================================================================
		// r = 3  stencil 1
		// =====================================================================
        
		IS_RHO(2)   = compute_third_order_smoothness_indicator(rho_coeff_31); 
		IS_U(2) = compute_third_order_smoothness_indicator(u_coeff_31);
		IS_V(2) = compute_third_order_smoothness_indicator(v_coeff_31);
		IS_P(2)     = compute_third_order_smoothness_indicator(p_coeff_31);
		IS_PHI(2)     = compute_third_order_smoothness_indicator(phi_coeff_31);
        
        
		// =====================================================================
		// r = 3  stencil 2
		// =====================================================================
        
		IS_RHO(3)   = compute_third_order_smoothness_indicator(rho_coeff_32); 
		IS_U(3) = compute_third_order_smoothness_indicator(u_coeff_32);
		IS_V(3) = compute_third_order_smoothness_indicator(v_coeff_32);
		IS_P(3)     = compute_third_order_smoothness_indicator(p_coeff_32);
   		IS_PHI(3)     = compute_third_order_smoothness_indicator(phi_coeff_32);
		// =====================================================================
		// r = 3  stencil 3
		// =====================================================================
        
		IS_RHO(4)   = compute_third_order_smoothness_indicator(rho_coeff_33); 
		IS_U(4) = compute_third_order_smoothness_indicator(u_coeff_33);
		IS_V(4) = compute_third_order_smoothness_indicator(v_coeff_33);
		IS_P(4)     = compute_third_order_smoothness_indicator(p_coeff_33);
		IS_PHI(4)     = compute_third_order_smoothness_indicator(phi_coeff_33);        
		// =====================================================================
		// r = 3  stencil 4
		// =====================================================================
        
		IS_RHO(5)   = compute_third_order_smoothness_indicator(rho_coeff_34); 
		IS_U(5) = compute_third_order_smoothness_indicator(u_coeff_34);
		IS_V(5) = compute_third_order_smoothness_indicator(v_coeff_34);
		IS_P(5)     = compute_third_order_smoothness_indicator(p_coeff_34);
		IS_PHI(5)     = compute_third_order_smoothness_indicator(phi_coeff_34);        
		// =====================================================================
		// r = 2  stencil 1
		// =====================================================================
		IS_RHO(6)   = compute_second_order_smoothness_indicator(rho_coeff_21); 
		IS_U(6) = compute_second_order_smoothness_indicator(u_coeff_21);
		IS_V(6) = compute_second_order_smoothness_indicator(v_coeff_21);
		IS_P(6)     = compute_second_order_smoothness_indicator(p_coeff_21);
		IS_PHI(6)     = compute_second_order_smoothness_indicator(phi_coeff_21);		
		// =====================================================================
		// r = 2  stencil 2
		// =====================================================================
		
		IS_RHO(7)   = compute_second_order_smoothness_indicator(rho_coeff_22); 
		IS_U(7) = compute_second_order_smoothness_indicator(u_coeff_22);
		IS_V(7) = compute_second_order_smoothness_indicator(v_coeff_22);
		IS_P(7)     = compute_second_order_smoothness_indicator(p_coeff_22);
		IS_PHI(7)     = compute_second_order_smoothness_indicator(phi_coeff_22);		
		// =====================================================================
		// r = 2  stencil 3
		// =====================================================================
		
		IS_RHO(8)   = compute_second_order_smoothness_indicator(rho_coeff_23); 
		IS_U(8) = compute_second_order_smoothness_indicator(u_coeff_23);
		IS_V(8) = compute_second_order_smoothness_indicator(v_coeff_23);
		IS_P(8)     = compute_second_order_smoothness_indicator(p_coeff_23);
		IS_PHI(8)     = compute_second_order_smoothness_indicator(phi_coeff_23);		
		// =====================================================================
		// r = 2  stencil 4
		// =====================================================================
		
		IS_RHO(9)   = compute_second_order_smoothness_indicator(rho_coeff_24); 
		IS_U(9) = compute_second_order_smoothness_indicator(u_coeff_24);
		IS_V(9) = compute_second_order_smoothness_indicator(v_coeff_24);
		IS_P(9)     = compute_second_order_smoothness_indicator(p_coeff_24);
		IS_PHI(9)     = compute_second_order_smoothness_indicator(phi_coeff_24);		
		// Combine the polynomials 
		sum_RHO = 0.0; sum_U = 0.0; sum_V = 0.0; sum_P = 0.0; sum_gamma = 0.0; 
		sum_PHI = 0.0;
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			
			w_RHO(j)   = gamma(j)/(std::pow((IS_RHO(j) + epsilon), p));
			w_U(j) = gamma(j)/(std::pow((IS_U(j) + epsilon), p));
			w_V(j) = gamma(j)/(std::pow((IS_V(j) + epsilon), p));
			w_P(j)     = gamma(j)/(std::pow((IS_P(j) + epsilon), p));
			w_PHI(j)     = gamma(j)/(std::pow((IS_PHI(j) + epsilon), p));
			
			sum_RHO += w_RHO(j); sum_U += w_U(j); sum_V += w_V(j); sum_P += w_P(j);
			sum_PHI += w_PHI(j);
			
			sum_gamma += gamma(j); 
		}
		
		// Normalize the weights 
		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			w_RHO(j) = w_RHO(j)/sum_RHO; 
			w_U(j) = w_U(j)/sum_U;
			w_V(j) = w_V(j)/sum_V;
			w_P(j) = w_P(j)/sum_P;
			w_PHI(j) = w_PHI(j)/sum_PHI;
			gamma(j) = gamma(j)/sum_gamma; 
		}
		// Density 
		
		coeffs_Rho[c](0) = rho0; 
		
		coeffs_Rho[c](1) = (w_RHO(0)/gamma(0))*(rho_coeff_4(0) - 
												
												gamma(1)*rho_coeff_3(0)  - 
												gamma(2)*rho_coeff_31(0) - 
												gamma(3)*rho_coeff_32(0) - 
												gamma(4)*rho_coeff_33(0) - 
												gamma(5)*rho_coeff_34(0) - 
												gamma(6)*rho_coeff_21(0) -
												gamma(7)*rho_coeff_22(0) -
												gamma(8)*rho_coeff_23(0) -
												gamma(9)*rho_coeff_24(0) ) + 
												
												w_RHO(1)*rho_coeff_3(0) +
												w_RHO(2)*rho_coeff_31(0) + 
												w_RHO(3)*rho_coeff_32(0) +
												w_RHO(4)*rho_coeff_33(0) +
												w_RHO(5)*rho_coeff_34(0) +
												w_RHO(6)*rho_coeff_21(0) +
												w_RHO(7)*rho_coeff_22(0) +
												w_RHO(8)*rho_coeff_23(0) +
												w_RHO(9)*rho_coeff_24(0) ;// u_x
		
		coeffs_Rho[c](2) = (w_RHO(0)/gamma(0))*(rho_coeff_4(1) - 
												
												gamma(1)*rho_coeff_3(1) - 
												gamma(2)*rho_coeff_31(1) - 
												gamma(3)*rho_coeff_32(1) -
												gamma(4)*rho_coeff_33(1) -
												gamma(5)*rho_coeff_34(1) - 
												gamma(6)*rho_coeff_21(1) -
												gamma(7)*rho_coeff_22(1) -
												gamma(8)*rho_coeff_23(1) -
												gamma(9)*rho_coeff_24(1) ) + 
												
												w_RHO(1)*rho_coeff_3(1) +
												w_RHO(2)*rho_coeff_31(1) +
												w_RHO(3)*rho_coeff_32(1) +
												w_RHO(4)*rho_coeff_33(1) +
												w_RHO(5)*rho_coeff_34(1) +
												w_RHO(6)*rho_coeff_21(1) +
												w_RHO(7)*rho_coeff_22(1) +
												w_RHO(8)*rho_coeff_23(1) +
												w_RHO(9)*rho_coeff_24(1) ; // u_x

		coeffs_Rho[c](3) = (w_RHO(0)/gamma(0))*(rho_coeff_4(2) - 
												
												gamma(1)*rho_coeff_3(2) - 
												gamma(2)*rho_coeff_31(2) - 
												gamma(3)*rho_coeff_32(2) - 
												gamma(4)*rho_coeff_33(2) - 
												gamma(5)*rho_coeff_34(2) ) + 
												
												w_RHO(1)*rho_coeff_3(2) +
												w_RHO(2)*rho_coeff_31(2) +
												w_RHO(3)*rho_coeff_32(2) + 
												w_RHO(4)*rho_coeff_33(2) + 
												w_RHO(5)*rho_coeff_34(2) ; // u_xx
												
		coeffs_Rho[c](4) = (w_RHO(0)/gamma(0))*(rho_coeff_4(3) - 
												
												gamma(1)*rho_coeff_3(3) - 
												gamma(2)*rho_coeff_31(3) - 
												gamma(3)*rho_coeff_32(3) - 
												gamma(4)*rho_coeff_33(3) - 
												gamma(5)*rho_coeff_34(3) ) + 
												
												w_RHO(1)*rho_coeff_3(3) +
												w_RHO(2)*rho_coeff_31(3) + 
												w_RHO(3)*rho_coeff_32(3) + 
												w_RHO(4)*rho_coeff_33(3) + 
												w_RHO(5)*rho_coeff_34(3) ; // u_yy
		
		coeffs_Rho[c](5) = (w_RHO(0)/gamma(0))*(rho_coeff_4(4) - 
												
												gamma(1)*rho_coeff_3(4) - 
												gamma(2)*rho_coeff_31(4) - 
												gamma(3)*rho_coeff_32(4) - 
												gamma(4)*rho_coeff_33(4) - 
												gamma(5)*rho_coeff_34(4) ) + 
												
												w_RHO(1)*rho_coeff_3(4) +
												w_RHO(2)*rho_coeff_31(4) + 
												w_RHO(3)*rho_coeff_32(4) + 
												w_RHO(4)*rho_coeff_33(4) + 
												w_RHO(5)*rho_coeff_34(4) ; // u_xy
		
		coeffs_Rho[c](6) = (w_RHO(0)/gamma(0))*(rho_coeff_4(5)); // u_xxx

		coeffs_Rho[c](7) = (w_RHO(0)/gamma(0))*(rho_coeff_4(6)); // u_yyy

		coeffs_Rho[c](8) = (w_RHO(0)/gamma(0))*(rho_coeff_4(7)); // u_xxy

		coeffs_Rho[c](9) = (w_RHO(0)/gamma(0))*(rho_coeff_4(8)); // u_xyy
		// x-momentum
		
		coeffs_U[c](0) = u0; 
		
		coeffs_U[c](1) = (w_U(0)/gamma(0))*(u_coeff_4(0) - 
												
												gamma(1)*u_coeff_3(0)  - 
												gamma(2)*u_coeff_31(0) - 
												gamma(3)*u_coeff_32(0) - 
												gamma(4)*u_coeff_33(0) - 
												gamma(5)*u_coeff_34(0) - 
												gamma(6)*u_coeff_21(0) -
												gamma(7)*u_coeff_22(0) -
												gamma(8)*u_coeff_23(0) -
												gamma(9)*u_coeff_24(0) ) + 
												
												w_U(1)*u_coeff_3(0) +
												w_U(2)*u_coeff_31(0) + 
												w_U(3)*u_coeff_32(0) +
												w_U(4)*u_coeff_33(0) +
												w_U(5)*u_coeff_34(0) +
												w_U(6)*u_coeff_21(0) +
												w_U(7)*u_coeff_22(0) +
												w_U(8)*u_coeff_23(0) +
												w_U(9)*u_coeff_24(0) ;// u_x
		
		coeffs_U[c](2) = (w_U(0)/gamma(0))*(u_coeff_4(1) - 
												
												gamma(1)*u_coeff_3(1) - 
												gamma(2)*u_coeff_31(1) - 
												gamma(3)*u_coeff_32(1) -
												gamma(4)*u_coeff_33(1) -
												gamma(5)*u_coeff_34(1) - 
												gamma(6)*u_coeff_21(1) -
												gamma(7)*u_coeff_22(1) -
												gamma(8)*u_coeff_23(1) -
												gamma(9)*u_coeff_24(1) ) + 
												
												w_U(1)*u_coeff_3(1) +
												w_U(2)*u_coeff_31(1) +
												w_U(3)*u_coeff_32(1) +
												w_U(4)*u_coeff_33(1) +
												w_U(5)*u_coeff_34(1) +
												w_U(6)*u_coeff_21(1) +
												w_U(7)*u_coeff_22(1) +
												w_U(8)*u_coeff_23(1) +
												w_U(9)*u_coeff_24(1) ; // u_x

		coeffs_U[c](3) = (w_U(0)/gamma(0))*(u_coeff_4(2) - 
												
												gamma(1)*u_coeff_3(2) - 
												gamma(2)*u_coeff_31(2) - 
												gamma(3)*u_coeff_32(2) - 
												gamma(4)*u_coeff_33(2) - 
												gamma(5)*u_coeff_34(2) ) + 
												
												w_U(1)*u_coeff_3(2) +
												w_U(2)*u_coeff_31(2) +
												w_U(3)*u_coeff_32(2) + 
												w_U(4)*u_coeff_33(2) + 
												w_U(5)*u_coeff_34(2) ; // u_xx
												
		coeffs_U[c](4) = (w_U(0)/gamma(0))*(u_coeff_4(3) - 
												
												gamma(1)*u_coeff_3(3) - 
												gamma(2)*u_coeff_31(3) - 
												gamma(3)*u_coeff_32(3) - 
												gamma(4)*u_coeff_33(3) - 
												gamma(5)*u_coeff_34(3) ) + 
												
												w_U(1)*u_coeff_3(3) +
												w_U(2)*u_coeff_31(3) + 
												w_U(3)*u_coeff_32(3) + 
												w_U(4)*u_coeff_33(3) + 
												w_U(5)*u_coeff_34(3) ; // u_yy
		
		coeffs_U[c](5) = (w_U(0)/gamma(0))*(u_coeff_4(4) - 
												
												gamma(1)*u_coeff_3(4) - 
												gamma(2)*u_coeff_31(4) - 
												gamma(3)*u_coeff_32(4) - 
												gamma(4)*u_coeff_33(4) - 
												gamma(5)*u_coeff_34(4) ) + 
												
												w_U(1)*u_coeff_3(4) +
												w_U(2)*u_coeff_31(4) + 
												w_U(3)*u_coeff_32(4) + 
												w_U(4)*u_coeff_33(4) + 
												w_U(5)*u_coeff_34(4) ; // u_xy
		
		coeffs_U[c](6) = (w_U(0)/gamma(0))*(u_coeff_4(5)); // u_xxx

		coeffs_U[c](7) = (w_U(0)/gamma(0))*(u_coeff_4(6)); // u_yyy

		coeffs_U[c](8) = (w_U(0)/gamma(0))*(u_coeff_4(7)); // u_xxy

		coeffs_U[c](9) = (w_U(0)/gamma(0))*(u_coeff_4(8)); // u_xyy
		// y-momentum
		
		coeffs_V[c](0) = v0; 
		
		coeffs_V[c](1) = (w_V(0)/gamma(0))*(v_coeff_4(0) - 
												
												gamma(1)*v_coeff_3(0)  - 
												gamma(2)*v_coeff_31(0) - 
												gamma(3)*v_coeff_32(0) - 
												gamma(4)*v_coeff_33(0) - 
												gamma(5)*v_coeff_34(0) - 
												gamma(6)*v_coeff_21(0) -
												gamma(7)*v_coeff_22(0) -
												gamma(8)*v_coeff_23(0) -
												gamma(9)*v_coeff_24(0) ) + 
												
												w_V(1)*v_coeff_3(0) +
												w_V(2)*v_coeff_31(0) + 
												w_V(3)*v_coeff_32(0) +
												w_V(4)*v_coeff_33(0) +
												w_V(5)*v_coeff_34(0) +
												w_V(6)*v_coeff_21(0) +
												w_V(7)*v_coeff_22(0) +
												w_V(8)*v_coeff_23(0) +
												w_V(9)*v_coeff_24(0) ;// v_x
		
		coeffs_V[c](2) = (w_V(0)/gamma(0))*(v_coeff_4(1) - 
												
												gamma(1)*v_coeff_3(1) - 
												gamma(2)*v_coeff_31(1) - 
												gamma(3)*v_coeff_32(1) -
												gamma(4)*v_coeff_33(1) -
												gamma(5)*v_coeff_34(1) - 
												gamma(6)*v_coeff_21(1) -
												gamma(7)*v_coeff_22(1) -
												gamma(8)*v_coeff_23(1) -
												gamma(9)*v_coeff_24(1) ) + 
												
												w_V(1)*v_coeff_3(1) +
												w_V(2)*v_coeff_31(1) +
												w_V(3)*v_coeff_32(1) +
												w_V(4)*v_coeff_33(1) +
												w_V(5)*v_coeff_34(1) +
												w_V(6)*v_coeff_21(1) +
												w_V(7)*v_coeff_22(1) +
												w_V(8)*v_coeff_23(1) +
												w_V(9)*v_coeff_24(1) ; // v_x

		coeffs_V[c](3) = (w_V(0)/gamma(0))*(v_coeff_4(2) - 
												
												gamma(1)*v_coeff_3(2) - 
												gamma(2)*v_coeff_31(2) - 
												gamma(3)*v_coeff_32(2) - 
												gamma(4)*v_coeff_33(2) - 
												gamma(5)*v_coeff_34(2) ) + 
												
												w_V(1)*v_coeff_3(2) +
												w_V(2)*v_coeff_31(2) +
												w_V(3)*v_coeff_32(2) + 
												w_V(4)*v_coeff_33(2) + 
												w_V(5)*v_coeff_34(2) ; // v_xx
												
		coeffs_V[c](4) = (w_V(0)/gamma(0))*(v_coeff_4(3) - 
												
												gamma(1)*v_coeff_3(3) - 
												gamma(2)*v_coeff_31(3) - 
												gamma(3)*v_coeff_32(3) - 
												gamma(4)*v_coeff_33(3) - 
												gamma(5)*v_coeff_34(3) ) + 
												
												w_V(1)*v_coeff_3(3) +
												w_V(2)*v_coeff_31(3) + 
												w_V(3)*v_coeff_32(3) + 
												w_V(4)*v_coeff_33(3) + 
												w_V(5)*v_coeff_34(3) ; // v_yy
		
		coeffs_V[c](5) = (w_V(0)/gamma(0))*(v_coeff_4(4) - 
												
												gamma(1)*v_coeff_3(4) - 
												gamma(2)*v_coeff_31(4) - 
												gamma(3)*v_coeff_32(4) - 
												gamma(4)*v_coeff_33(4) - 
												gamma(5)*v_coeff_34(4) ) + 
												
												w_V(1)*v_coeff_3(4) +
												w_V(2)*v_coeff_31(4) + 
												w_V(3)*v_coeff_32(4) + 
												w_V(4)*v_coeff_33(4) + 
												w_V(5)*v_coeff_34(4) ; // v_xy
		
		coeffs_V[c](6) = (w_V(0)/gamma(0))*(v_coeff_4(5)); // v_xxx

		coeffs_V[c](7) = (w_V(0)/gamma(0))*(v_coeff_4(6)); // v_yyy

		coeffs_V[c](8) = (w_V(0)/gamma(0))*(v_coeff_4(7)); // v_xxy

		coeffs_V[c](9) = (w_V(0)/gamma(0))*(v_coeff_4(8)); // v_xyy
		// Total energy  
		
		coeffs_P[c](0) = p0; 
		
		coeffs_P[c](1) = (w_P(0)/gamma(0))*(p_coeff_4(0) - 
												
												gamma(1)*p_coeff_3(0)  - 
												gamma(2)*p_coeff_31(0) - 
												gamma(3)*p_coeff_32(0) - 
												gamma(4)*p_coeff_33(0) - 
												gamma(5)*p_coeff_34(0) - 
												gamma(6)*p_coeff_21(0) -
												gamma(7)*p_coeff_22(0) -
												gamma(8)*p_coeff_23(0) -
												gamma(9)*p_coeff_24(0) ) + 
												
												w_P(1)*p_coeff_3(0) +
												w_P(2)*p_coeff_31(0) + 
												w_P(3)*p_coeff_32(0) +
												w_P(4)*p_coeff_33(0) +
												w_P(5)*p_coeff_34(0) +
												w_P(6)*p_coeff_21(0) +
												w_P(7)*p_coeff_22(0) +
												w_P(8)*p_coeff_23(0) +
												w_P(9)*p_coeff_24(0) ;// v_x
		
		coeffs_P[c](2) = (w_P(0)/gamma(0))*(p_coeff_4(1) - 
												
												gamma(1)*p_coeff_3(1) - 
												gamma(2)*p_coeff_31(1) - 
												gamma(3)*p_coeff_32(1) -
												gamma(4)*p_coeff_33(1) -
												gamma(5)*p_coeff_34(1) - 
												gamma(6)*p_coeff_21(1) -
												gamma(7)*p_coeff_22(1) -
												gamma(8)*p_coeff_23(1) -
												gamma(9)*p_coeff_24(1) ) + 
												
												w_P(1)*p_coeff_3(1) +
												w_P(2)*p_coeff_31(1) +
												w_P(3)*p_coeff_32(1) +
												w_P(4)*p_coeff_33(1) +
												w_P(5)*p_coeff_34(1) +
												w_P(6)*p_coeff_21(1) +
												w_P(7)*p_coeff_22(1) +
												w_P(8)*p_coeff_23(1) +
												w_P(9)*p_coeff_24(1) ; // v_x

		coeffs_P[c](3) = (w_P(0)/gamma(0))*(p_coeff_4(2) - 
												
												gamma(1)*p_coeff_3(2) - 
												gamma(2)*p_coeff_31(2) - 
												gamma(3)*p_coeff_32(2) - 
												gamma(4)*p_coeff_33(2) - 
												gamma(5)*p_coeff_34(2) ) + 
												
												w_P(1)*p_coeff_3(2) +
												w_P(2)*p_coeff_31(2) +
												w_P(3)*p_coeff_32(2) + 
												w_P(4)*p_coeff_33(2) + 
												w_P(5)*p_coeff_34(2) ; // v_xx
												
		coeffs_P[c](4) = (w_P(0)/gamma(0))*(p_coeff_4(3) - 
												
												gamma(1)*p_coeff_3(3) - 
												gamma(2)*p_coeff_31(3) - 
												gamma(3)*p_coeff_32(3) - 
												gamma(4)*p_coeff_33(3) - 
												gamma(5)*p_coeff_34(3) ) + 
												
												w_P(1)*p_coeff_3(3) +
												w_P(2)*p_coeff_31(3) + 
												w_P(3)*p_coeff_32(3) + 
												w_P(4)*p_coeff_33(3) + 
												w_P(5)*p_coeff_34(3) ; // v_yy
		
		coeffs_P[c](5) = (w_P(0)/gamma(0))*(p_coeff_4(4) - 
												
												gamma(1)*p_coeff_3(4) - 
												gamma(2)*p_coeff_31(4) - 
												gamma(3)*p_coeff_32(4) - 
												gamma(4)*p_coeff_33(4) - 
												gamma(5)*p_coeff_34(4) ) + 
												
												w_P(1)*p_coeff_3(4) +
												w_P(2)*p_coeff_31(4) + 
												w_P(3)*p_coeff_32(4) + 
												w_P(4)*p_coeff_33(4) + 
												w_P(5)*p_coeff_34(4) ; // v_xy
		
		coeffs_P[c](6) = (w_P(0)/gamma(0))*(p_coeff_4(5)); // u_xxx

		coeffs_P[c](7) = (w_P(0)/gamma(0))*(p_coeff_4(6)); // u_yyy

		coeffs_P[c](8) = (w_P(0)/gamma(0))*(p_coeff_4(7)); // u_xxy

		coeffs_P[c](9) = (w_P(0)/gamma(0))*(p_coeff_4(8)); // u_xyy
		// Phi 
		coeffs_Phi[c](0) = phi0; 
		
		coeffs_Phi[c](1) = (w_PHI(0)/gamma(0))*(phi_coeff_4(0) - 
												
												gamma(1)*phi_coeff_3(0)  - 
												gamma(2)*phi_coeff_31(0) - 
												gamma(3)*phi_coeff_32(0) - 
												gamma(4)*phi_coeff_33(0) - 
												gamma(5)*phi_coeff_34(0) - 
												gamma(6)*phi_coeff_21(0) -
												gamma(7)*phi_coeff_22(0) -
												gamma(8)*phi_coeff_23(0) -
												gamma(9)*phi_coeff_24(0) ) + 
												
												w_PHI(1)*phi_coeff_3(0) +
												w_PHI(2)*phi_coeff_31(0) + 
												w_PHI(3)*phi_coeff_32(0) +
												w_PHI(4)*phi_coeff_33(0) +
												w_PHI(5)*phi_coeff_34(0) +
												w_PHI(6)*phi_coeff_21(0) +
												w_PHI(7)*phi_coeff_22(0) +
												w_PHI(8)*phi_coeff_23(0) +
												w_PHI(9)*phi_coeff_24(0) ;// v_x
		
		coeffs_Phi[c](2) = (w_PHI(0)/gamma(0))*(phi_coeff_4(1) - 
												
												gamma(1)*phi_coeff_3(1) - 
												gamma(2)*phi_coeff_31(1) - 
												gamma(3)*phi_coeff_32(1) -
												gamma(4)*phi_coeff_33(1) -
												gamma(5)*phi_coeff_34(1) - 
												gamma(6)*phi_coeff_21(1) -
												gamma(7)*phi_coeff_22(1) -
												gamma(8)*phi_coeff_23(1) -
												gamma(9)*phi_coeff_24(1) ) + 
												
												w_PHI(1)*phi_coeff_3(1) +
												w_PHI(2)*phi_coeff_31(1) +
												w_PHI(3)*phi_coeff_32(1) +
												w_PHI(4)*phi_coeff_33(1) +
												w_PHI(5)*phi_coeff_34(1) +
												w_PHI(6)*phi_coeff_21(1) +
												w_PHI(7)*phi_coeff_22(1) +
												w_PHI(8)*phi_coeff_23(1) +
												w_PHI(9)*phi_coeff_24(1) ; // v_x

		coeffs_Phi[c](3) = (w_PHI(0)/gamma(0))*(phi_coeff_4(2) - 
												
												gamma(1)*phi_coeff_3(2) - 
												gamma(2)*phi_coeff_31(2) - 
												gamma(3)*phi_coeff_32(2) - 
												gamma(4)*phi_coeff_33(2) - 
												gamma(5)*phi_coeff_34(2) ) + 
												
												w_PHI(1)*phi_coeff_3(2) +
												w_PHI(2)*phi_coeff_31(2) +
												w_PHI(3)*phi_coeff_32(2) + 
												w_PHI(4)*phi_coeff_33(2) + 
												w_PHI(5)*phi_coeff_34(2) ; // v_xx
												
		coeffs_Phi[c](4) = (w_PHI(0)/gamma(0))*(phi_coeff_4(3) - 
												
												gamma(1)*phi_coeff_3(3) - 
												gamma(2)*phi_coeff_31(3) - 
												gamma(3)*phi_coeff_32(3) - 
												gamma(4)*phi_coeff_33(3) - 
												gamma(5)*phi_coeff_34(3) ) + 
												
												w_PHI(1)*phi_coeff_3(3) +
												w_PHI(2)*phi_coeff_31(3) + 
												w_PHI(3)*phi_coeff_32(3) + 
												w_PHI(4)*phi_coeff_33(3) + 
												w_PHI(5)*phi_coeff_34(3) ; // v_yy
		
		coeffs_Phi[c](5) = (w_PHI(0)/gamma(0))*(phi_coeff_4(4) - 
												
												gamma(1)*phi_coeff_3(4) - 
												gamma(2)*phi_coeff_31(4) - 
												gamma(3)*phi_coeff_32(4) - 
												gamma(4)*phi_coeff_33(4) - 
												gamma(5)*phi_coeff_34(4) ) + 
												
												w_PHI(1)*phi_coeff_3(4) +
												w_PHI(2)*phi_coeff_31(4) + 
												w_PHI(3)*phi_coeff_32(4) + 
												w_PHI(4)*phi_coeff_33(4) + 
												w_PHI(5)*phi_coeff_34(4) ; // v_xy
		
		coeffs_Phi[c](6) = (w_PHI(0)/gamma(0))*(phi_coeff_4(5)); // u_xxx

		coeffs_Phi[c](7) = (w_PHI(0)/gamma(0))*(phi_coeff_4(6)); // u_yyy

		coeffs_Phi[c](8) = (w_PHI(0)/gamma(0))*(phi_coeff_4(7)); // u_xxy

		coeffs_Phi[c](9) = (w_PHI(0)/gamma(0))*(phi_coeff_4(8)); // u_xyy

		// check bounds 
//        if(!Use_ader)

        for (unsigned int f = 0; f < faces_per_cell; ++f) {

            face_quadrature_point_1 = Cell[c].face_quadrature_point1(f);
            face_quadrature_point_2 = Cell[c].face_quadrature_point2(f);

            W1(0) = evaluate_weno_polynomial(coeffs_Rho[c], WENO_poly_consts[c], face_quadrature_point_1, h);
            W1(1) = evaluate_weno_polynomial(coeffs_U[c], WENO_poly_consts[c], face_quadrature_point_1, h);
            W1(2) = evaluate_weno_polynomial(coeffs_V[c], WENO_poly_consts[c], face_quadrature_point_1, h);
            W1(3) = evaluate_weno_polynomial(coeffs_P[c], WENO_poly_consts[c], face_quadrature_point_1, h);
            W1(4) = evaluate_weno_polynomial(coeffs_Phi[c], WENO_poly_consts[c], face_quadrature_point_1, h);

            W2(0) = evaluate_weno_polynomial(coeffs_Rho[c], WENO_poly_consts[c], face_quadrature_point_2, h);
            W2(1) = evaluate_weno_polynomial(coeffs_U[c], WENO_poly_consts[c], face_quadrature_point_2, h);
            W2(2) = evaluate_weno_polynomial(coeffs_V[c], WENO_poly_consts[c], face_quadrature_point_2, h);
            W2(3) = evaluate_weno_polynomial(coeffs_P[c], WENO_poly_consts[c], face_quadrature_point_2, h);
            W2(4) = evaluate_weno_polynomial(coeffs_Phi[c], WENO_poly_consts[c], face_quadrature_point_2, h);

            if (W1(0) < 0.0 || W1(3) < 0.0 || W1(4) < 0.0 || W1(4) > 1.0 || W2(0) < 0.0 || W2(3) < 0.0 || W2(4) < 0.0 || W2(4) > 1.0) {

                coeffs_Rho[c](0) = rho0;  coeffs_U[c](0) = u0;    coeffs_V[c](0) = v0;    coeffs_P[c](0) = p0; 
				coeffs_Phi[c](0) = phi0; 

                for (unsigned int i = 1; i < 10; i++) {
                  coeffs_Rho[c](i) = 0.0; 
                  coeffs_U[c](i) = 0.0;
                  coeffs_V[c](i) = 0.0;
                  coeffs_P[c](i) = 0.0;
                  coeffs_Phi[c](i) = 0.0;
                }

            }
        }
		//std::cout<<"c: "<<c<<std::endl;
    } // End of cell loop 
//	MPI_Barrier(MPI_COMM_WORLD);

} // End of function 
