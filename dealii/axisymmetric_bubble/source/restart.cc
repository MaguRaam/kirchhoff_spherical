#include "../include/Weno432.h"

void Weno4_2D::restart() {

	if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
		std::ofstream statusR_S("../restart/Status.dat", std::ios::out);
	    statusR_S.flags(std::ios::dec |std::ios::scientific);
		statusR_S.precision(16);

	    if ( !(statusR_S.is_open()) ) {
	        std::cerr << "Unable to open restart folder" << std::endl; 
	        std::exit(EXIT_FAILURE);
	    }		

    	statusR_S << time <<std::endl;
    	statusR_S << int(0.0) << std::endl;
		statusR_S.close();

	}
   
    std::ofstream fout_RHO("../restart/RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 5) + ".dat");
    if ( !(fout_RHO.is_open()) ) {
        std::cerr << "Unable to open restart folder" << std::endl; 
        std::exit(EXIT_FAILURE);
    }
    
    fout_RHO.flags( std::ios::dec | std::ios::scientific );
	fout_RHO.precision(16);

    fout_RHO << time << std::endl; 

	unsigned int g_i;

	for (unsigned int c = 0; c < n_locally_cells; ++c) {

		g_i = local_to_global_index_map[c];

	        fout_RHO << local_RHO(g_i ) << std::endl; 
	        fout_RHO << local_RHO_U(g_i ) << std::endl;
	        fout_RHO << local_RHO_V(g_i ) << std::endl;
	        fout_RHO << local_E(g_i ) << std::endl;
	        fout_RHO << local_PHI(g_i ) << std::endl;
	}
    
    fout_RHO.close(); 

	if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
		std::ofstream statusR_E("../restart/Status.dat",std::ios::out);
	    statusR_E.flags(std::ios::dec |std::ios::scientific);
		statusR_E.precision(16);

	    if ( !(statusR_E.is_open()) ) {
	        std::cerr << "Unable to open restart folder" << std::endl; 
	        std::exit(EXIT_FAILURE);
	    }		

    	statusR_E << time <<std::endl;
    	statusR_E << int(1.0) << std::endl;
		statusR_E.close();

	}
   
}

void Weno4_2D::restart_r() {

	if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
		std::ofstream statusR_S("../restart_r/Status.dat",std::ios::out);
	    statusR_S.flags(std::ios::dec |std::ios::scientific);
		statusR_S.precision(16);

	    if ( !(statusR_S.is_open()) ) {
	        std::cerr << "Unable to open restart_r folder" << std::endl; 
	        std::exit(EXIT_FAILURE);
	    }		

    	statusR_S << time <<std::endl;
    	statusR_S << int(0.0) <<std::endl;
		statusR_S.close();

	}
	unsigned int g_i;   
    std::ofstream fout_RHO("../restart_r/RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 5) + ".dat");
    if ( !(fout_RHO.is_open()) ) {
        std::cerr << "Unable to open restart_r folder" << std::endl; 
        std::exit(EXIT_FAILURE);
    }
    
    fout_RHO.flags( std::ios::dec | std::ios::scientific );
	fout_RHO.precision(16);

    fout_RHO << time << std::endl; 

	for (unsigned int c = 0; c < n_locally_cells; ++c) {

		g_i = local_to_global_index_map[c];

        fout_RHO << local_RHO(g_i ) << std::endl; 
        fout_RHO << local_RHO_U(g_i ) << std::endl;
        fout_RHO << local_RHO_V(g_i ) << std::endl;
        fout_RHO << local_E(g_i ) << std::endl;
        fout_RHO << local_PHI(g_i ) << std::endl;
	}
    
    fout_RHO.close(); 

	if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) {
		std::ofstream statusR_E("../restart_r/Status.dat",std::ios::out);
	    statusR_E.flags(std::ios::dec |std::ios::scientific);
		statusR_E.precision(16);

	    if ( !(statusR_E.is_open()) ) {
	        std::cerr << "Unable to open restart_r folder" << std::endl; 
	        std::exit(EXIT_FAILURE);
	    }		

    	statusR_E << time <<std::endl;
    	statusR_E << int(1.0) <<std::endl;
		statusR_E.close();

	}
    
}

