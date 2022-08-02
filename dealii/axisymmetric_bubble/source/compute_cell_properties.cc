#include "../include/Weno432.h" 

void Weno4_2D::compute_cell_properties() {

    // Iterate over all the cells 
//    auto start = std::chrono::system_clock::now();
    
    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
/*
    // For obtaining cell quadrature points and jxws
    QGauss<2> quadrature_formula(2);
    FEValues<2> fv_values (mapping, fv, quadrature_formula, update_quadrature_points | update_JxW_values);
    
    // For obtaining face normal vectors
	QGauss<2-1> face_quadrature_formula(2);
	FEFaceValues<2> fv_face_values (mapping, fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors | update_JxW_values);

    // For obtaining face normal vectors at face center
	QGauss<2-1> face_quadrature_formula_center(1);
	FEFaceValues<2> fv_face_values_center (mapping, fv, face_quadrature_formula_center, update_quadrature_points | update_normal_vectors);
*/

    // For obtaining cell quadrature points and jxws
    QGauss<2> quadrature_formula(2);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);
    
    // For obtaining face normal vectors
	QGauss<2-1> face_quadrature_formula(2);
	FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors | update_JxW_values);

    // For obtaining face normal vectors at face center
	QGauss<2-1> face_quadrature_formula_center(1);
	FEFaceValues<2> fv_face_values_center (fv, face_quadrature_formula_center, update_quadrature_points | update_normal_vectors);


    Point<2> face_quadrature_point_1;     
    Point<2> face_quadrature_point_2;
    
    double S_f; // Volume of the cell and surface area of the face
	
	double nx1, ny1;   // Face normal vectors
	double nx2, ny2; 
	
	Tensor<1,2> face_normal_vector1; // Face normal vector
	Tensor<1,2> face_normal_vector2; // Face normal vector
	Tensor<1,2> face_center_normal_vector; // Face normal vector
    
    double area;
    Point<2> cell_quadrature_points[4];
    double jxws[4];     
    Point<2> face_quadrature_points1[4];
    Point<2> face_quadrature_points2[4];
    double face_jxws_1[4];      
    double face_jxws_2[4];      
    double face_normal_x1[4]; 
    double face_normal_x2[4]; 
    double face_normal_y1[4]; 
    double face_normal_y2[4];
    Point<2> face_center_quadrature_points[4];
    double face_center_normal_x[4]; 
    double face_center_normal_y[4];
    double surface_area[4];

	double h_old,h_min_local; 
	h_min_local = 1e6; 
    
	for (; cell != endc; ++cell) {

        cell->get_dof_indices(local_dof_indices);

		if (is_ghost_cell[local_dof_indices[0] ] ){
           
	        area = 0.0; 

			fv_values.reinit(cell);
			for (unsigned int i=0; i<fv_values.n_quadrature_points; ++i){
				area += fv_values.JxW (i);
	            cell_quadrature_points[i] = fv_values.quadrature_point(i);
	            jxws[i] =  fv_values.JxW (i);            
	        }

			h_old = std::sqrt(area);
	
			if (h_old < h_min_local) {
				h_min_local = h_old; 
			}

        
	        for (unsigned int f = 0; f < 4; f++) {
                
				fv_face_values.reinit(cell, f);
    	        fv_face_values_center.reinit(cell, f);                
    	            
    	        face_normal_vector1 = fv_face_values.normal_vector(0);
    	        nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
    	            
    	        face_normal_vector2 = fv_face_values.normal_vector(1);
    	        nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
    	            
    	        face_quadrature_point_1 = fv_face_values.quadrature_point(0); 
    	        face_quadrature_point_2 = fv_face_values.quadrature_point(1); 
	
    	        face_center_normal_vector = fv_face_values_center.normal_vector(0);  
  	
    	        S_f = 0.0; 
    	   		for (unsigned int i=0; i<fv_face_values.n_quadrature_points; ++i)		
		            S_f += fv_face_values.JxW (i);
    	            
   	            // Fill the values 
    	            
    	        face_quadrature_points1[f] = face_quadrature_point_1;  
    	        face_quadrature_points2[f] = face_quadrature_point_2;  
    	        face_normal_x1[f] = nx1;
    	        face_normal_x2[f] = nx2; 
    	        face_normal_y1[f] = ny1; 
    	        face_normal_y2[f] = ny2; 
				face_jxws_1[f] =  fv_face_values.JxW (0);
				face_jxws_2[f] =  fv_face_values.JxW (1);
    	        face_center_quadrature_points[f] = fv_face_values_center.quadrature_point(0); 
    	        face_center_normal_x[f] = face_center_normal_vector[0]; 
    	        face_center_normal_y[f] = face_center_normal_vector[1]; 
    	        surface_area[f] = S_f;
	
			}
				
			unsigned int c = global_to_local_index_map[local_dof_indices[0] ];
            Cell[c].reinit(area, cell_quadrature_points, jxws, face_quadrature_points1, face_quadrature_points2,
face_jxws_1, face_jxws_2, face_normal_x1, face_normal_x2, face_normal_y1,
 face_normal_y2,face_center_quadrature_points, face_center_normal_x, face_center_normal_y, surface_area); 
		}
	}

	h_min = Utilities::MPI::min (h_min_local, MPI_COMM_WORLD);
    pcout<<"h_min: "<<h_min<<std::endl;
	pcout<<"compute cell properties"<<std::endl;
/*
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "timer.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Time taken to compute cell properties = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}
*/
}
