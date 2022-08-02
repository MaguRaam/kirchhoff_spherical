#include "../include/Weno432.h"


//  Evaluate time step using the CFL condition 

void Weno4_2D::compute_time_step_based_on_cfl_number() {

	Vector<double> u(n_locally_cells);
	Vector<double> v(n_locally_cells);
//	Vector<double> a(n_locally_cells);
    
    Vector<double> U_(5); Vector<double> W(5); 

	unsigned int g_i;
	double a;
	double gamma_, p_inf_;
    
	for (unsigned int c = 0; c < n_locally_cells; ++c) {

		g_i = local_to_global_index_map[c];        
        U_(0) = RHO(g_i); U_(1) = RHO_U(g_i); U_(2) = RHO_V(g_i); U_(3) = E(g_i); U_(4) = PHI(g_i); 

		gamma_ = claw.gamma(U_(4));
		p_inf_ = claw.P_inf(U_(4), gamma_);
        
        claw.conserved_to_primitive(U_, gamma_, p_inf_, W); 

        a = std::sqrt(gamma_ * ( W(3) + p_inf_) /W(0));

        u(c) = a + std::abs(W(1));
        v(c) = a + std::abs(W(2));

	}
	
	double u_max_local = u.linfty_norm();
	double v_max_local = v.linfty_norm();

	double u_max = Utilities::MPI::max (u_max_local, MPI_COMM_WORLD);
	double v_max = Utilities::MPI::max (v_max_local, MPI_COMM_WORLD);

	dt = (cfl*h_min)/(sqrt(2.0)*sqrt(u_max*u_max + v_max*v_max));

	if((time + dt)>finalTime) {
		dt = finalTime- time;
	}
} 
