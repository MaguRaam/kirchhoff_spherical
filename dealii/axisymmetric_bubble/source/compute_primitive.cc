#include "../include/Weno432.h"


//  Evaluate time step using the CFL condition 

void Weno4_2D::compute_primitive() {
    
    Vector<double> u(5); Vector<double> W(5); 
	double gamma_, p_inf_;
	unsigned int g_i;

//	if(Utilities::MPI::this_mpi_process(mpi_communicator) == 0)     
	for (unsigned int c = 0; c < n_store_cell; ++c) {

		g_i = local_to_global_index_map[c];        
        u(0) = RHO(g_i); u(1) = RHO_U(g_i); u(2) = RHO_V(g_i); u(3) = E(g_i); u(4) = PHI(g_i); 

		gamma_ = claw.gamma(u(4));
		p_inf_ = claw.P_inf(u(4), gamma_);
        
        claw.conserved_to_primitive(u, gamma_, p_inf_, W);

		Rho(c) = W(0);
		U(c) = W(1);
		V(c) = W(2);
		P(c) = W(3);
		Phi(c) = W(4);
//		std::cout<<"g_i: "<<g_i<<"\tc: "<<c<<std::endl;
	}
//	MPI_Barrier(mpi_communicator);
//	std::cout<<"compute_primitive: \n";

/*
	if(Utilities::MPI::this_mpi_process(mpi_communicator) == 1) 
    
	for (unsigned int c = 0; c < n_store_cell; ++c) {

		g_i = local_to_global_index_map[c];        
        U(0) = RHO(g_i); U(1) = RHO_U(g_i); U(2) = RHO_V(g_i); U(3) = E(g_i); U(4) = PHI(g_i); 

		gamma_ = claw.gamma(U(4));
		p_inf_ = claw.P_inf(U(4), gamma_);
        
        claw.conserved_to_primitive(U, gamma_, p_inf_, W);

		Rho(c) = W(0);
		U(c) = W(1);
		V(c) = W(2);
		P(c) = W(3);
		Phi(c) = W(4);

//		std::cout<<"g_i: "<<g_i<<"\tc: "<<c<<std::endl;
	}
*/
//	std::cout<<"compute_primitive: \n"<<Utilities::MPI::this_mpi_process(mpi_communicator)<<std::endl;
//	MPI_Barrier(mpi_communicator);


} 
