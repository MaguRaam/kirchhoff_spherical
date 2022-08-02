#include "../include/Riemann_Solvers.h"
#include "../include/Headers.h"
#include "../include/Exceptions.h"


// Function definitions

// 1. primitive_to_conserved
Vector<double> primitive_to_conserved(Vector<double> W) {
	Vector<double> U(4);

	double e = W[3]/((GAMMA - 1.0)*W[0]);   // Internal energy, e = p/((gamma-1)*rho)
	double k = 0.5*(W[1]*W[1] + W[2]*W[2]); // Kinetic energy, k = (1/2)*(u^2 + v^2)

	U[0] = W[0];
	U[1] = W[0]*W[1];
	U[2] = W[0]*W[2];
	U[3] = W[0]*(k + e);

	return U;
}

// 2. conserved_to_primitive
Vector<double> conserved_to_primitive(Vector<double> U) {

	Vector<double> W(4);

	W[0] = U[0];
	W[1] = U[1]/U[0];
	W[2] = U[2]/U[0];
	W[3] = (GAMMA-1.0)*(U[3] - 0.5*U[0]*(W[1]*W[1]+ W[2]*W[2]));

	return W;
}


//3. compute_flux_from_primitive_variable
Vector<double> compute_flux_from_primitive_variable(Vector<double> W, double nx, double ny) {

		Vector<double> FV(4);
		Vector<double> Q(4);

		 Q = primitive_to_conserved(W);

		double un; // normal velocity

		un = (Q[1]*nx + Q[2]*ny) ;
      FV[0] = Q[0]*un ;
      FV[1] = Q[0]*Q[1]*un + Q[3]*nx ;
      FV[2] = Q[0]*Q[2]*un + Q[3]*ny ;
      FV[3] = ( GAMMA*(Q[3])/(GAMMA - 1.0) + 0.5*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2]) )*un ;

      return FV;
}


// Lax Friedrichs Riemann Solver

void Usual_To_Conservative(double*& Var_Q, double* Var_U, double gamma_l, double P_Infty_l) {
        Var_Q[0] = Var_U[0] ; // rho
        Var_Q[1] = Var_U[0]*Var_U[1] ; // rho*u
        Var_Q[2] = Var_U[0]*Var_U[2] ; // rho*v
        Var_Q[3] = Var_U[3]/(gamma_l - 1.0) + gamma_l*P_Infty_l/(gamma_l - 1.0) + 0.5*Var_U[0]*(Var_U[1]*Var_U[1] + Var_U[2]*Var_U[2]) ; // E
}

void FluxFunction(double*& FV, double* Q, double gamma_l, double P_Infty_l, double nx, double ny) {
	double un; // normal velocity
	un = (Q[1]*nx + Q[2]*ny) ;
        FV[0] = Q[0]*un ;
        FV[1] = Q[0]*Q[1]*un + Q[3]*nx ;
        FV[2] = Q[0]*Q[2]*un + Q[3]*ny ;
        FV[3] = ( gamma_l*(Q[3] + P_Infty_l)/(gamma_l - 1.0) + 0.5*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2]) )*un ;
}

// Input is Primitive Vector
Vector<double> local_Lax_Friedrichs_riemann_solver(Vector <double> UL, Vector <double> UR, double nx, double ny, Point<2> P, bool boundary) {
    
    if (UL[0] < 0.0 || UR[0] < 0.0 || UL[3] < 0.0 || UR[3] < 0.0) {
        ThrowNegativePressureDensityError(UL, UR, P, boundary); 
    }

	double U_L[] = {UL[0], UL[1], UL[2], UL[3]};
   double U_R[] = {UR[0], UR[1], UR[2], UR[3]};
   Vector<double> F(4);
	double flux[4];
 	double S_star;
 	double P_Infty_L = 0.0;
 	double P_Infty_R = 0.0;
 	int i ;
	double rho_L, rho_R, u_L, u_R, v_L, v_R, P_L, P_R, c_L, c_R ;
	double un_L, un_R, S_L, S_R ;
	//double ut_L, ut_R;
	double *FL, *FR, *U ;

	U = new double[4] ; FL = new double[4] ; FR = new double[4] ;

	rho_L = U_L[0] ; u_L = U_L[1] ; v_L = U_L[2] ; P_L = U_L[3] ;
  rho_R = U_R[0] ; u_R = U_R[1] ; v_R = U_R[2] ; P_R = U_R[3] ;
	un_L = u_L*nx + v_L*ny ; //ut_L = -u_L*ny + v_L*nx ;
	un_R = u_R*nx + v_R*ny ; //ut_R = -u_R*ny + v_R*nx ;

  c_L = sqrt(GAMMA*(P_L + P_Infty_L)/rho_L) ; c_R = sqrt(GAMMA*(P_R + P_Infty_R)/rho_R) ;

	
	S_L = std::max(fabs(un_R - c_R), fabs(un_L - c_L)) ;
	S_R = std::max(fabs(un_L + c_L), fabs(un_R + c_R)) ;
	S_star = std::max(fabs(S_L),fabs(S_R)) ;

	// Now compute the left right and starred fluxes.
	Usual_To_Conservative(U,U_L,GAMMA,P_Infty_L) ;
	FluxFunction(FL,U_L,GAMMA,P_Infty_L,nx,ny) ;
        for(i = 0 ; i < 4 ; i++) flux[i] = 0.5*FL[i] + 0.5*S_star*U[i] ;

	Usual_To_Conservative(U,U_R,GAMMA,P_Infty_R) ;
	FluxFunction(FR,U_R,GAMMA,P_Infty_R,nx,ny);
        for(i = 0 ; i < 4 ; i++) flux[i] += 0.5*FR[i] - 0.5*S_star*U[i] ;

    F[0] = flux[0]; F[1] = flux[1]; F[2] = flux[2]; F[3] = flux[3];

	delete [] U ; delete [] FL ; delete [] FR ;
	return F;
}


// Return HLLC flux at midpoint 0.0 ;
Vector<double> HLLC_riemann_solver(Vector<double> UL, Vector<double> UR, double nx, double ny, Point<2> P, bool boundary) {
    
    if (UL[0] < 0.0 || UR[0] < 0.0 || UL[3] < 0.0 || UR[3] < 0.0) {
        ThrowNegativePressureDensityError(UL, UR, P, boundary); 
    }
    
    
	Vector<double> F(4);
    double flux[4];
    int i ;
    double U_L[] = {UL[0], UL[1], UL[2], UL[3]};
    double U_R[] = {UR[0], UR[1], UR[2], UR[3]};
	double rho_L, rho_R, u_L, u_R, v_L, v_R, P_L, P_R, c_L, c_R, E_L, E_R ;
	double un_L, un_R, ut_L, ut_R ;
	double un, ut ;
	double S_L, S_R ;
	double *UL_star, *UR_star, *FL, *FR, *FL_star, *FR_star, *U ;
    double P_Infty_ = 0.0;
    double S_star; 
	
	UL_star = new double[4] ; UR_star = new double[4] ; U = new double[4] ;
	FL = new double[4] ; FR = new double[4] ; FL_star = new double[4] ; FR_star = new double[4] ;

	rho_L = U_L[0] ; u_L = U_L[1] ; v_L = U_L[2] ; P_L = U_L[3] ;
        rho_R = U_R[0] ; u_R = U_R[1] ; v_R = U_R[2] ; P_R = U_R[3] ;
	un_L = u_L*nx + v_L*ny ; ut_L = -u_L*ny + v_L*nx ;
	un_R = u_R*nx + v_R*ny ; ut_R = -u_R*ny + v_R*nx ;

    	c_L = sqrt(GAMMA*(P_L + P_Infty_)/rho_L) ; c_R = sqrt(GAMMA*(P_R + P_Infty_)/rho_R) ;


	E_L = (P_L + GAMMA*P_Infty_)/(GAMMA - 1.0) + 0.5*rho_L*u_L*u_L + 0.5*rho_L*v_L*v_L ;
	E_R = (P_R + GAMMA*P_Infty_)/(GAMMA - 1.0) + 0.5*rho_R*u_R*u_R + 0.5*rho_R*v_R*v_R ;

	S_L = std::min((un_R - c_R), (un_L - c_L)) ; 
	S_R = std::max((un_L + c_L), (un_R + c_R)) ;
	S_star = ( P_R - P_L + rho_L*un_L*(S_L-un_L) - rho_R*un_R*(S_R - un_R) )/(rho_L*(S_L - un_L) - rho_R*(S_R - un_R)) ;
	//P_Star = P_L + rho_L*(S_star - un_L)*(S_L - un_L) ;

	// Now compute the left right and starred fluxes for HLLC.
	FluxFunction(FL,U_L,GAMMA,P_Infty_,nx,ny) ; FluxFunction(FR,U_R,GAMMA,P_Infty_,nx,ny);
	
	UL_star[0] = rho_L*(S_L - un_L)/(S_L - S_star) ;  
	un = S_star ; ut = ut_L ;
	UL_star[1] = UL_star[0]*(un*nx - ut*ny) ;  
	UL_star[2] = UL_star[0]*(un*ny + ut*nx) ;
	UL_star[3] = UL_star[0]*( (E_L/rho_L) + (S_star - un_L)*(S_star + P_L/(rho_L*(S_L - un_L)) ) )	;  
	
	UR_star[0] = rho_R*(S_R - un_R)/(S_R - S_star) ;
	un = S_star ; ut = ut_R ;
	UR_star[1] = UR_star[0]*(un*nx - ut*ny) ;
	UR_star[2] = UR_star[0]*(un*ny + ut*nx) ;
	UR_star[3] = UR_star[0]*( (E_R/rho_R) + (S_star - un_R)*(S_star + P_R/(rho_R*(S_R - un_R)) ) ) ;
	
	Usual_To_Conservative(U,U_L,GAMMA,P_Infty_) ;
	U[4] = U_L[4] ; 
	for(i = 0 ; i < 4 ; i++) FL_star[i] = FL[i] + S_L*(UL_star[i] - U[i]) ; 

	Usual_To_Conservative(U,U_R,GAMMA,P_Infty_) ;
	U[4] = U_R[4] ; 
	for(i = 0 ; i < 4 ; i++) FR_star[i] = FR[i] + S_R*(UR_star[i] - U[i]) ; 

	if( S_L > 0.0 ) {
		for(i = 0 ; i < 4 ; i++) flux[i] = FL[i] ; 
	} else if((S_star >= 0.0) && (S_L < 0.0)) {
		for(i = 0 ; i < 4 ; i++) flux[i] = FL_star[i] ; 
	} else if((S_star < 0.0) && (S_R >= 0.0)) {
		for(i = 0 ; i < 4 ; i++) flux[i] = FR_star[i] ; 
	} else if(S_R < 0.0) {
		for(i = 0 ; i < 4 ; i++) flux[i] = FR[i] ; 
	}

	delete [] UL_star ; delete [] UR_star ; delete [] FL ; delete [] FR ; delete [] FL_star ; delete [] FR_star ;
	delete [] U ;

    F[0] = flux[0]; F[1] = flux[1]; F[2] = flux[2]; F[3] = flux[3];

    return F;
}


Vector<double> rotated_HLLC_riemann_solver(Vector<double> UL, Vector<double> UR, double nx, double ny, Point<2> P, bool boundary) {
 
    Vector<double> flux1(4);
    Vector<double> flux2(4);
    double flux[4];
    Vector<double> F(4);
    double U_L[] = {UL[0], UL[1], UL[2], UL[3]};
    double U_R[] = {UR[0], UR[1], UR[2], UR[3]};

    int i ;
	double alpha1, alpha2, n1x, n1y, n2x, n2y, u_L, u_R, v_L, v_R, du, dv, dq ;

	u_L = U_L[1] ; v_L = U_L[2] ; 
    u_R = U_R[1] ; v_R = U_R[2] ; 
	du = u_R - u_L ; dv = v_R - v_L ; 
	dq = sqrt(du*du + dv*dv) ;
	
	if(dq < 1.0E-10) { n1x = nx ; n1y = ny ; }
	else { n1x = du/dq ; n1y = dv/dq ; }

	alpha1 = (n1x*nx + n1y*ny) ;
	if(alpha1 < 0) { n1x = -n1x ; n1y = -n1y ; alpha1 = -alpha1 ; }
	n2x = -n1y ; n2y = n1x ;
	alpha2 = (n2x*nx + n2y*ny) ; 
	if(alpha2 < 0) { n2x = -n2x ; n2y = -n2y ; alpha2 = -alpha2 ; }

    flux1 = HLLC_riemann_solver(UL, UR, n1x, n1y, P, boundary);
	flux2 = HLLC_riemann_solver(UL, UR, n2x, n2y, P, boundary);

	for(i = 0 ; i < 4 ; i++) flux[i] = alpha1*flux1[i] + alpha2*flux2[i] ;
	
    F[0] = flux[0]; F[1] = flux[1]; F[2] = flux[2]; F[3] = flux[3];

    return F;
}

