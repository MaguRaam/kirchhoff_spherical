#include "../include/Weno432.h"

// Put everything together - solve the actual problem

void Weno4_2D::run()
{

	auto start = std::chrono::system_clock::now();

	pcout << "=============================================================="
		  << std::endl
		  << "Running with "
		  << "PETSc"
		  << " on "
		  << Utilities::MPI::n_mpi_processes(mpi_communicator)
		  << " MPI rank(s)..." << std::endl;

	make_grid();

	triangulation.refine_global(n_refinement);

	std::size_t memory = triangulation.memory_consumption();

	pcout << "memory in GB: " << memory / std::pow(2, 30) << std::endl;

	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
		std::ofstream fout_convergence;
		fout_convergence.flags(std::ios::dec | std::ios::scientific);
		fout_convergence.precision(7);

		const std::string filename = "log.dat";
		fout_convergence.open(filename, std::ios::in | std::ios::out | std::ios::app);

		fout_convergence << "Number of active cells: " << triangulation.n_active_cells() << std::endl
						 << "memory in GB: " << memory / std::pow(2, 30) << std::endl;
		fout_convergence.close();
	}

	pcout << "Number of active cells: "
		  << triangulation.n_active_cells()
		  << std::endl
		  << "memory in GB: " << memory / std::pow(2, 30) << std::endl;
	pcout << "============================" << std::endl;

	auto start_grid_output = std::chrono::system_clock::now();

	if(!RESTART) output_grid();

	auto start_setup_system = std::chrono::system_clock::now();

	setup_system();

	get_quad_points(); //for kirchhoff solver

	auto start_allocate_memory = std::chrono::system_clock::now();

	allocate_memory();

	auto start_compute_cell_properties = std::chrono::system_clock::now();

	compute_cell_properties();

	auto start_compute_weno_polynomial_constants = std::chrono::system_clock::now();

	compute_weno_polynomial_constants();

	auto start_initialize = std::chrono::system_clock::now();

	initialize();

	auto start_precompute_matrices = std::chrono::system_clock::now();

	precompute_matrices();

	auto start_compute_IS_constants = std::chrono::system_clock::now();

	compute_IS_constants();

	auto start_solve_ader = std::chrono::system_clock::now();

	std::chrono::duration<double> grid = start_grid_output - start;
	std::chrono::duration<double> grid_output = start_setup_system - start_grid_output;
	std::chrono::duration<double> setup_system = start_allocate_memory - start_setup_system;
	std::chrono::duration<double> allocate_memory = start_compute_cell_properties - start_allocate_memory;
	std::chrono::duration<double> compute_cell_properties = start_compute_weno_polynomial_constants - start_compute_cell_properties;
	std::chrono::duration<double> compute_weno_polynomial_constants = start_initialize - start_compute_weno_polynomial_constants;
	std::chrono::duration<double> initialize = start_precompute_matrices - start_initialize;
	std::chrono::duration<double> precompute_matrices = start_compute_IS_constants - start_precompute_matrices;
	std::chrono::duration<double> compute_IS_constants = start_solve_ader - start_compute_IS_constants;
	std::chrono::duration<double> pre_ader = start_solve_ader - start;

	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
		std::ofstream fout_convergence;
		fout_convergence.flags(std::ios::dec | std::ios::scientific);
		fout_convergence.precision(7);

		const std::string filename = "timer.dat";
		fout_convergence.open(filename, std::ios::in | std::ios::out | std::ios::app);

		fout_convergence << "total cell: " << triangulation.n_active_cells() << std::endl
						 << "cell per core: " << n_locally_cells << std::endl
						 << "relevant cell: " << n_relevant_cells << std::endl
						 << "grid generation : " << grid.count() << std::endl
						 << "grid_output: " << grid_output.count() << std::endl
						 << "setup_system: " << setup_system.count() << std::endl
						 << "allocate_memory: " << allocate_memory.count() << std::endl
						 << "compute_cell_properties: " << compute_cell_properties.count() << std::endl
						 << "compute_weno_polynomial_constants: " << compute_weno_polynomial_constants.count() << std::endl
						 << "initialize: " << initialize.count() << std::endl
						 << "precompute_matrices: " << precompute_matrices.count() << std::endl
						 << "compute_IS_constants: " << compute_IS_constants.count() << std::endl
						 << "pre_ader: " << pre_ader.count() << std::endl;
		fout_convergence.close();
	}

	start_solve_ader = std::chrono::system_clock::now();

	//	solve_ader();
		solve_ssprk33();

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> solve_ader = end - start_solve_ader;
	std::chrono::duration<double> total = end - start;

	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
		std::ofstream fout_convergence;
		fout_convergence.flags(std::ios::dec | std::ios::scientific);
		fout_convergence.precision(7);

		const std::string filename = "timer.dat";
		fout_convergence.open(filename, std::ios::in | std::ios::out | std::ios::app);

		fout_convergence << "total cell: " << triangulation.n_active_cells() << std::endl
						 << "cell per core: " << n_locally_cells << std::endl
						 << "relevant cell: " << n_relevant_cells << std::endl
						 << "solve_ader: " << solve_ader.count() << std::endl
						 << "Total time = " << total.count() << std::endl;
		fout_convergence.close();
	}

	computing_timer.print_summary();
}
