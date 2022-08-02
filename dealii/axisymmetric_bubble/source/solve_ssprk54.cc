#include "../include/Weno432.h"


// Update the solution using SSPRK(5,4) method

void Weno4_2D::solve_ssprk54() {

	Use_ader = false;

    // Using SSPRK(5,4)

    auto start = std::chrono::system_clock::now();

    unsigned int count = 0;

    Vector<double> RHO_old(n_locally_cells);
    Vector<double> RHO_U_old(n_locally_cells);
	Vector<double> RHO_V_old(n_locally_cells);
	Vector<double> E_old(n_locally_cells); 
	Vector<double> PHI_old(n_locally_cells); 
    
    Vector<double> RHO_2(n_locally_cells);
    Vector<double> RHO_U_2(n_locally_cells);
	Vector<double> RHO_V_2(n_locally_cells);
	Vector<double> E_2(n_locally_cells); 
	Vector<double> PHI_2(n_locally_cells); 
    
    Vector<double> RHO_3(n_locally_cells);
    Vector<double> RHO_U_3(n_locally_cells);
	Vector<double> RHO_V_3(n_locally_cells);
	Vector<double> E_3(n_locally_cells);
	Vector<double> PHI_3(n_locally_cells); 
    
    Vector<double> rhs1_3(n_locally_cells);
    Vector<double> rhs2_3(n_locally_cells);
	Vector<double> rhs3_3(n_locally_cells);
	Vector<double> rhs4_3(n_locally_cells); 
	Vector<double> rhs5_3(n_locally_cells); 

	compute_time_step_based_on_cfl_number();

	unsigned int output_count;
	output_count = 0.2 / dt ;

	unsigned int g_i ;

	while (time < finalTime) {

	    auto start_ssprk54 = std::chrono::system_clock::now();

		compute_time_step_based_on_cfl_number();

		if (count%output_count == 0 ) {
            output_results();
        }

		if (count%100 == 0 ) {
            restart();
        }

		if ((count+5)%100 == 0 ) {
            restart_r();
        }

		time += dt;
        
    	if ( (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) && ( ( count%50 == 0 ) || std::fabs(time - finalTime) < 1e-8) ){
		    std::ofstream fout_convergence ;
    	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
    	    fout_convergence.precision(7) ;
	
    	    const std::string filename = "log.dat";
		    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    		fout_convergence <<"time = " << time <<",\t dt: "<<dt<< ",\t Final time = " << finalTime<< std::endl;
    	    fout_convergence.close();
		}

		pcout<<"time = " << time <<",\t dt: "<<dt<< ",\t Final time = " << finalTime<< std::endl;

        count++;

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

			RHO_old(c) = RHO(g_i);
			RHO_U_old(c) = RHO_U(g_i);
			RHO_V_old(c) = RHO_V(g_i);
			E_old(c) = E(g_i);	
			PHI_old(c) = PHI(g_i);	
		}
        
        compute_rhs();

		// SSPRK Stage 1
		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

            local_RHO(g_i) = local_RHO(g_i) + 0.391752226571890*dt*rhs1(c);
            local_RHO_U(g_i) = local_RHO_U(g_i) + 0.391752226571890*dt*rhs2(c);
            local_RHO_V(g_i) = local_RHO_V(g_i) + 0.391752226571890*dt*rhs3(c);
            local_E(g_i) = local_E(g_i) + 0.391752226571890*dt*rhs4(c);
            local_PHI(g_i) = local_PHI(g_i) + 0.391752226571890*dt*rhs5(c);

        }

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);
		local_PHI.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;
		PHI = local_PHI;

		// SSPRK Stage 2		
        
		compute_rhs();

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

            local_RHO(g_i)   = 0.444370493651235*RHO_old(c)   + 0.555629506348765*local_RHO(g_i)   + 0.368410593050371*dt*rhs1(c);
            local_RHO_U(g_i) = 0.444370493651235*RHO_U_old(c) + 0.555629506348765*local_RHO_U(g_i) + 0.368410593050371*dt*rhs2(c);
            local_RHO_V(g_i) = 0.444370493651235*RHO_V_old(c) + 0.555629506348765*local_RHO_V(g_i) + 0.368410593050371*dt*rhs3(c);
            local_E(g_i)     = 0.444370493651235*E_old(c)     + 0.555629506348765*local_E(g_i)     + 0.368410593050371*dt*rhs4(c);
            local_PHI(g_i)     = 0.444370493651235*PHI_old(c)     + 0.555629506348765*local_PHI(g_i)     + 0.368410593050371*dt*rhs5(c);
        }

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);
		local_PHI.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;
		PHI = local_PHI;

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

			RHO_2(c) = RHO(g_i);
			RHO_U_2(c) = RHO_U(g_i);
			RHO_V_2(c) = RHO_V(g_i);
			E_2(c) = E(g_i);
			PHI_2(c) = PHI(g_i);
		}

		// SSPRK Stage 3

		compute_rhs();

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

            local_RHO(g_i)   = 0.620101851488403*RHO_old(c)   + 0.379898148511597*local_RHO(g_i)   + 0.251891774271694*dt*rhs1(c);
            local_RHO_U(g_i) = 0.620101851488403*RHO_U_old(c) + 0.379898148511597*local_RHO_U(g_i) + 0.251891774271694*dt*rhs2(c);
            local_RHO_V(g_i) = 0.620101851488403*RHO_V_old(c) + 0.379898148511597*local_RHO_V(g_i) + 0.251891774271694*dt*rhs3(c);
            local_E(g_i)     = 0.620101851488403*E_old(c)     + 0.379898148511597*local_E(g_i)     + 0.251891774271694*dt*rhs4(c);
            local_PHI(g_i)     = 0.620101851488403*PHI_old(c)     + 0.379898148511597*local_PHI(g_i)     + 0.251891774271694*dt*rhs5(c);
		}

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);
		local_PHI.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;
		PHI = local_PHI;
		
		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

			RHO_3(c) = RHO(g_i);
			RHO_U_3(c) = RHO_U(g_i);
			RHO_V_3(c) = RHO_V(g_i);
			E_3(c) = E(g_i);
			PHI_3(c) = PHI(g_i);
		}
        
		
        // SSPRK Stage 4

		compute_rhs();
        
        rhs1_3 = rhs1; rhs2_3 = rhs2; rhs3_3 = rhs3; rhs4_3 = rhs4; rhs5_3 = rhs5;

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

            local_RHO(g_i)   = 0.178079954393132*RHO_old(c)   + 0.821920045606868*local_RHO(g_i)   + 0.544974750228521*dt*rhs1(c);
            local_RHO_U(g_i) = 0.178079954393132*RHO_U_old(c) + 0.821920045606868*local_RHO_U(g_i) + 0.544974750228521*dt*rhs2(c);
            local_RHO_V(g_i) = 0.178079954393132*RHO_V_old(c) + 0.821920045606868*local_RHO_V(g_i) + 0.544974750228521*dt*rhs3(c);
            local_E(g_i)     = 0.178079954393132*E_old(c)     + 0.821920045606868*local_E(g_i)     + 0.544974750228521*dt*rhs4(c);
            local_PHI(g_i)     = 0.178079954393132*PHI_old(c)     + 0.821920045606868*local_PHI(g_i)     + 0.544974750228521*dt*rhs5(c);
		}

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);
		local_PHI.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;
		PHI = local_PHI;
		
		//  SSPRK Stage 5
		
        compute_rhs();
		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

            local_RHO(g_i)   = 0.517231671970585*RHO_2(c)   + 0.096059710526147*RHO_3(c)   + 0.386708617503269*local_RHO(g_i) +  
                        0.226007483236906*dt*rhs1(c)+ 0.063692468666290*dt*rhs1_3(c);
            
            local_RHO_U(g_i) = 0.517231671970585*RHO_U_2(c) + 0.096059710526147*RHO_U_3(c) + 0.386708617503269*local_RHO_U(g_i) +
                        0.226007483236906*dt*rhs2(c)+ 0.063692468666290*dt*rhs2_3(c);
            
            local_RHO_V(g_i) = 0.517231671970585*RHO_V_2(c) + 0.096059710526147*RHO_V_3(c) + 0.386708617503269*local_RHO_V(g_i) +
                        0.226007483236906*dt*rhs3(c)+ 0.063692468666290*dt*rhs3_3(c);
            
            local_E(g_i)     = 0.517231671970585*E_2(c)     + 0.096059710526147*E_3(c)     + 0.386708617503269*local_E(g_i) + 
                        0.226007483236906*dt*rhs4(c)+ 0.063692468666290*dt*rhs4_3(c);

            local_PHI(g_i)     = 0.517231671970585*PHI_2(c)     + 0.096059710526147*PHI_3(c)     + 0.386708617503269*local_PHI(g_i) + 
                        0.226007483236906*dt*rhs5(c)+ 0.063692468666290*dt*rhs5_3(c);

        }

	    auto start_com = std::chrono::system_clock::now();

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);
		local_PHI.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;
		PHI = local_PHI;

	    auto end_com = std::chrono::system_clock::now();
	    std::chrono::duration<double> elapsed_seconds_com = end_com - start_com;
	    std::chrono::duration<double> elapsed_seconds_ssprk54 = end_com - start_ssprk54;

    	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0 && count == 2){
		    std::ofstream fout_convergence ;
    	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
    	    fout_convergence.precision(7) ;
	
    	    const std::string filename = "timer.dat";
		    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);
	
    		fout_convergence << "time taken by data transfer of ssprk54 = " << elapsed_seconds_com.count() << std::endl;
    		fout_convergence << "time taken by 1 step of ssprk54 = " << elapsed_seconds_ssprk54.count() << std::endl;
    	    fout_convergence.close();
		}
        
    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "timer.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "time taken by ssprk54 = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}
    
    output_results();
    restart(); 
} 
