#include "../include/Weno432.h"

// Solve using SSPRK(3,3)



void Weno4_2D::solve_ssprk33()
{
	pcout << "solve by ssprk33: " << std::endl;

	//	Use_ader = false;

	auto start = std::chrono::system_clock::now();

	unsigned int count = 0;

	Vector<double> RHO_old(n_locally_cells);
	Vector<double> RHO_U_old(n_locally_cells);
	Vector<double> RHO_V_old(n_locally_cells);
	Vector<double> E_old(n_locally_cells);
	Vector<double> PHI_old(n_locally_cells);

	compute_primitive();
	reconstruct();
	compute_time_step_based_on_cfl_number();
	unsigned int output_count;
	output_count = 0.0001 / dt;
	//	if( output_count == 0 ) output_count = 10;
	pcout << "reconstruct: " << std::endl;

	unsigned int g_i;

	// open pressure and time file:
	std::ofstream pfile, tfile, pofile;

	if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
	{
		pfile.open("../kirchhoff/data/pressure.dat");
		tfile.open("../kirchhoff/data/time.dat");
		pfile.flags(std::ios::dec | std::ios::scientific);
		tfile.flags(std::ios::dec | std::ios::scientific);
		pfile.precision(16);
		tfile.precision(16);
	}

	while (time < finalTime)
	{

		// write pressure on semicircle arc (Kirchhoff surface):
		if (count % 10 == 0)
		{

			get_pressure(); // compute pressure on arc:

			// write data from 0th process:
			if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
			{

				// write time data:
				pfile << "t = " << time << "\n";
				tfile << time << "\n";

				for (unsigned int i = 0; i < q_global.size(); ++i)
					pfile << q_global[i].p << "\t" << q_global[i].pr << "\n";

			}

			//write pressure at observer point:
			auto cell = GridTools::find_active_cell_around_point(mapping, dof_handler, po).first;

			if (cell != dof_handler.end() && cell->is_locally_owned()){

				cell->get_dof_indices(local_dof_indices);
    			unsigned int local_index = global_to_local_index_map[local_dof_indices[0]];
    		    double observer_pressure = evaluate_weno_polynomial(coeffs_P[local_index], WENO_poly_consts[local_index], po, Cell[local_index].h()) - p_inf;

    		    pofile.open("../kirchhoff/p_weno.dat", std::ios_base::app);
    		    pofile.flags(std::ios::dec | std::ios::scientific);
    		    pofile.precision(16);

    		    pofile << time << "\t" << observer_pressure << "\n";
			}


		}

		auto start_ssprk33 = std::chrono::system_clock::now();

		compute_time_step_based_on_cfl_number();
		/*
				if (count < 100) {
					dt = dt/(-0.09*static_cast<double>(count) + 10.0);
				}
		*/

		if (count % 200 == 0)
		{
			output_results();
		}

		if (count % 300 == 0)
		{
			restart();
		}

		if ((count + 5) % 600 == 0)
		{
			restart_r();
		}

		time += dt;

		if ((Utilities::MPI::this_mpi_process(mpi_communicator) == 0) && ((count % 50 == 0) || std::fabs(time - finalTime) < 1e-8))
		{
			std::ofstream fout_convergence;
			fout_convergence.flags(std::ios::dec | std::ios::scientific);
			fout_convergence.precision(7);

			const std::string filename = "log.dat";
			fout_convergence.open(filename, std::ios::in | std::ios::out | std::ios::app);

			fout_convergence << "time = " << time << ",\t dt: " << dt << ",\t Final time = " << finalTime << std::endl;
			fout_convergence.close();
		}

		pcout << "time = " << time << ",\t dt: " << dt << ",\t Final time = " << finalTime << std::endl;

		for (unsigned int c = 0; c < n_locally_cells; ++c)
		{
			g_i = local_to_global_index_map[c];
			RHO_old(c) = RHO(g_i);
			RHO_U_old(c) = RHO_U(g_i);
			RHO_V_old(c) = RHO_V(g_i);
			E_old(c) = E(g_i);
			PHI_old(c) = PHI(g_i);
		}

		// SSPRK Step 1

		compute_rhs();

		//		pcout<<"1st step"<<std::endl;

		for (unsigned int c = 0; c < n_locally_cells; ++c)
		{

			g_i = local_to_global_index_map[c];

			local_RHO(g_i) = local_RHO(g_i) + dt * rhs1(c);
			local_RHO_U(g_i) = local_RHO_U(g_i) + dt * rhs2(c);
			local_RHO_V(g_i) = local_RHO_V(g_i) + dt * rhs3(c);
			local_E(g_i) = local_E(g_i) + dt * rhs4(c);
			local_PHI(g_i) = local_PHI(g_i) + dt * rhs5(c);
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

		// SSPRK Step 2

		compute_rhs();

		for (unsigned int c = 0; c < n_locally_cells; ++c)
		{

			g_i = local_to_global_index_map[c];

			local_RHO(g_i) = (3. / 4.) * RHO_old(c) + (1. / 4.) * local_RHO(g_i) + (1. / 4.) * dt * rhs1(c);
			local_RHO_U(g_i) = (3. / 4.) * RHO_U_old(c) + (1. / 4.) * local_RHO_U(g_i) + (1. / 4.) * dt * rhs2(c);
			local_RHO_V(g_i) = (3. / 4.) * RHO_V_old(c) + (1. / 4.) * local_RHO_V(g_i) + (1. / 4.) * dt * rhs3(c);
			local_E(g_i) = (3. / 4.) * E_old(c) + (1. / 4.) * local_E(g_i) + (1. / 4.) * dt * rhs4(c);
			local_PHI(g_i) = (3. / 4.) * PHI_old(c) + (1. / 4.) * local_PHI(g_i) + (1. / 4.) * dt * rhs5(c);
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

		// SSPRK Step 3

		compute_rhs();

		for (unsigned int c = 0; c < n_locally_cells; ++c)
		{

			g_i = local_to_global_index_map[c];

			local_RHO(g_i) = (1. / 3.) * RHO_old(c) + (2. / 3.) * local_RHO(g_i) + (2. / 3.) * dt * rhs1(c);
			local_RHO_U(g_i) = (1. / 3.) * RHO_U_old(c) + (2. / 3.) * local_RHO_U(g_i) + (2. / 3.) * dt * rhs2(c);
			local_RHO_V(g_i) = (1. / 3.) * RHO_V_old(c) + (2. / 3.) * local_RHO_V(g_i) + (2. / 3.) * dt * rhs3(c);
			local_E(g_i) = (1. / 3.) * E_old(c) + (2. / 3.) * local_E(g_i) + (2. / 3.) * dt * rhs4(c);
			local_PHI(g_i) = (1. / 3.) * PHI_old(c) + (2. / 3.) * local_PHI(g_i) + (2. / 3.) * dt * rhs5(c);
		}

		auto start_com = std::chrono::system_clock::now();

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);
		local_PHI.compress(VectorOperation::insert);

		auto end_com = std::chrono::system_clock::now();

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;
		PHI = local_PHI;

		auto end_trans = std::chrono::system_clock::now();

		std::chrono::duration<double> elapsed_seconds_com = end_com - start_com;
		std::chrono::duration<double> elapsed_seconds_trans = end_trans - end_com;
		std::chrono::duration<double> elapsed_seconds_ssprk33 = end_trans - start_ssprk33;

		if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0 && (count < 4))
		{
			std::ofstream fout_convergence;
			fout_convergence.flags(std::ios::dec | std::ios::scientific);
			fout_convergence.precision(7);

			const std::string filename = "timer.dat";
			fout_convergence.open(filename, std::ios::in | std::ios::out | std::ios::app);

			fout_convergence << "time taken by data compression of ssprk33 = " << elapsed_seconds_com.count() << std::endl;
			fout_convergence << "time taken by data transfer of ssprk33 = " << elapsed_seconds_trans.count() << std::endl;
			fout_convergence << "time taken by 1 step of ssprk33 = " << elapsed_seconds_ssprk33.count() << std::endl;
			fout_convergence.close();
		}
		count++;
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
		std::ofstream fout_convergence;
		fout_convergence.flags(std::ios::dec | std::ios::scientific);
		fout_convergence.precision(7);

		const std::string filename = "timer.dat";
		fout_convergence.open(filename, std::ios::in | std::ios::out | std::ios::app);

		fout_convergence << "time taken by ssprk33 = " << elapsed_seconds.count() << std::endl;
		fout_convergence.close();
	}

	output_results();
	restart();
}
