/*
 * claw.cc
 *      Author: sunder
 */

#include "../include/claw.h"

void ErrNegativePressureDensity(double d, double p, double x, double y) {

	if (d < 0.0 || p < 0.0) {
	
		std::cerr << "Negative Pressure/Density" << std::endl;
        std::cerr << "Density = " << d << ", Pressure = " << p << std::endl;
        std::cerr << "at (" << x << ", " << y << ")" << std::endl;
		
		std::exit(1); 
	}
}

// Conservation law definitions

//----------------------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------------------


Conservation_Law::Conservation_Law() :
	prandtl_no(0.0),
	GAMMA_L(0.0),
	GAMMA_G(0.0),
	MU_L(0.0),
	MU_G(0.0),
	P_INF_L(0.0),
	P_INF_G(0.0)
{}

//----------------------------------------------------------------------------
// Constructor taking specific heat at constant pressure, specific heat at constant
// volume and Prandtl number as inputs
//----------------------------------------------------------------------------

Conservation_Law::Conservation_Law(double pr, double gamma_l, double gamma_g, double mu_l, double mu_g, double p_inf_l, double p_inf_g) :
	prandtl_no(pr),
	GAMMA_L(gamma_l),
	GAMMA_G(gamma_g),
	MU_L(mu_l),
	MU_G(mu_g),
	P_INF_L(p_inf_l),
	P_INF_G(p_inf_g)
{
//	gamma() = cp/cv;
}

//----------------------------------------------------------------------------
// Get ratio of specific heats
//----------------------------------------------------------------------------

double Conservation_Law::gamma() const {
	return 1.4;
}

// phi = 0 indicates air phase---> gamma_Fluid = GAMMA_G

double Conservation_Law::gamma(const Vector<double>& U) {
	double gamma_Fluid = 1.0 + (GAMMA_L-1.0)*(GAMMA_G-1.0)/(GAMMA_G*U(4) + GAMMA_L*(1.0 - U(4)) - 1.0 ) ;
//	if (gamma_Fluid < GAMMA_G) gamma_Fluid = GAMMA_G;
//	if (gamma_Fluid > GAMMA_L) gamma_Fluid = GAMMA_L;
	return gamma_Fluid;
}

double Conservation_Law::gamma(const double phi) {
	double gamma_Fluid = 1.0 + (GAMMA_L-1.0)*(GAMMA_G-1.0)/(GAMMA_G*phi + GAMMA_L*(1.0 - phi) - 1.0 ) ;
//	if (gamma_Fluid < GAMMA_G) gamma_Fluid = GAMMA_G;
//	if (gamma_Fluid > GAMMA_L) gamma_Fluid = GAMMA_L;
	return gamma_Fluid;
}

//----------------------------------------------------------------------------
// Get P_inf of mixture
//----------------------------------------------------------------------------


double Conservation_Law::P_inf(const Vector<double>& U) {
	double gamma_Fluid = 1.0 + (GAMMA_L-1.0)*(GAMMA_G-1.0)/(GAMMA_G*U(4) + GAMMA_L*(1.0 - U(4)) - 1.0 ) ;
//	if (gamma_Fluid < GAMMA_G) gamma_Fluid = GAMMA_G;
//	if (gamma_Fluid > GAMMA_L) gamma_Fluid = GAMMA_L;
	double P_Infty_Fluid = (gamma_Fluid-1.0)*( P_INF_L*U(4)*GAMMA_L/(GAMMA_L-1.0) + P_INF_G*(1.0 - U(4))*GAMMA_G/(GAMMA_G-1.0) )/gamma_Fluid ;
//	if (P_Infty_Fluid < P_INF_G) P_Infty_Fluid = P_INF_G;
//	if (P_Infty_Fluid > P_INF_L) P_Infty_Fluid = P_INF_L;
	return P_Infty_Fluid;
}


double Conservation_Law::P_inf(const Vector<double>& U, const double gamma_Fluid) {
	double P_Infty_Fluid = (gamma_Fluid-1.0)*( P_INF_L*U(4)*GAMMA_L/(GAMMA_L-1.0) + P_INF_G*(1.0 - U(4))*GAMMA_G/(GAMMA_G-1.0) )/gamma_Fluid ;
//	if (P_Infty_Fluid < P_INF_G) P_Infty_Fluid = P_INF_G;
//	if (P_Infty_Fluid > P_INF_L) P_Infty_Fluid = P_INF_L;	
	return P_Infty_Fluid;
}

double Conservation_Law::P_inf(const double phi, const double gamma_Fluid) {
	double P_Infty_Fluid = (gamma_Fluid-1.0)*( P_INF_L*phi*GAMMA_L/(GAMMA_L-1.0) + P_INF_G*(1.0 - phi)*GAMMA_G/(GAMMA_G-1.0) )/gamma_Fluid ;
//	if (P_Infty_Fluid < P_INF_G) P_Infty_Fluid = P_INF_G;
//	if (P_Infty_Fluid > P_INF_L) P_Infty_Fluid = P_INF_L;
	return P_Infty_Fluid;
}

//----------------------------------------------------------------------------
// Find the viscosity of the gas at the given state U
//----------------------------------------------------------------------------

double Conservation_Law::viscosity(const Vector<double>& U)  {

	// Use Sutherland's law
	/*
	double T = Temperature(U);
	static const double mu_0 = 0.1;
	static const double T_0 = 1.0;
	static const double beta =
	double T = Tempe	static const double beta =rature(U);
	double T = Temperature(U);
	double T = Temperature(U); 1.5;
	static const double s = 1.0;

	double mu = mu_0*std::pow(T/T_0, beta)*((T_0 + s)/(T + s));
	*/
	double mu = MU_L * U(4) + MU_G * (1.0 - U(4));

	return mu;
}

//----------------------------------------------------------------------------
// Get Prandtl number
//----------------------------------------------------------------------------

double Conservation_Law::Pr() const {
	return prandtl_no;
}

//----------------------------------------------------------------------------
// Find heat conductivity of the gas at the given state U
//----------------------------------------------------------------------------

double Conservation_Law::heat_conductivity(const Vector<double>& U) {

	return ((gamma()/(gamma() - 1.0))*viscosity(U)/prandtl_no);
}

//----------------------------------------------------------------------------
// Given a conserved variable, find the pressure of the gas
//----------------------------------------------------------------------------

double Conservation_Law::Pressure(const Vector<double>& U) {

	return ( (gamma()-1.0)*(U[3] - 0.5*((U[1]*U[1] + U[2]*U[2])/U[0])) );
}

//----------------------------------------------------------------------------
// Given a conserved variable, find temperature of the gas
//----------------------------------------------------------------------------

double Conservation_Law::Temperature(const Vector<double>& U) {


	return ( Pressure(U)/(U[0]) );
}

//----------------------------------------------------------------------------
// Given a conserved variable, find the speed of sound
//----------------------------------------------------------------------------

double Conservation_Law::speed_of_sound(const Vector<double>& U) {

	return std::sqrt(gamma()*Pressure(U)/U[0]);
}

//----------------------------------------------------------------------------
// Get conserved variable from a primitive variable
//----------------------------------------------------------------------------

void Conservation_Law::primitive_to_conserved(const Vector<double>& W, Vector<double>& U) const {
/*
  Assert(W.size() == 5,
         ExcDimensionMismatch(W.size(), 5) );
  Assert(U.size() == 5,
         ExcDimensionMismatch(U.size(), 5) );
*/
	//double p = W[3]*W[0]*R(); // p = rho*R*T

    double e = W[3]/((gamma() - 1.0)*W[0]);
    double k = 0.5*(W[1]*W[1] + W[2]*W[2]);

    U[0] = W[0];
    U[1] = W[0]*W[1];
    U[2] = W[0]*W[2];
    U[3] = W[0]*(k + e);
    U[4] = W[4];
}

void Conservation_Law::primitive_to_conserved(const Vector<double>& W, double gamma, double p_inf, Vector<double>& U) const {
/*
  Assert(W.size() == 5,
         ExcDimensionMismatch(W.size(), 5) );
  Assert(U.size() == 5,
         ExcDimensionMismatch(U.size(), 5) );
*/
	//double p = W[3]*W[0]*R(); // p = rho*R*T

    double e = W[3]/(gamma - 1.0) + gamma*p_inf/(gamma - 1.0);
    double k = W[0]*0.5*(W[1]*W[1] + W[2]*W[2]);
    U[0] = W[0];
    U[1] = W[0]*W[1];
    U[2] = W[0]*W[2];
    U[3] = k + e;
    U[4] = W[4];
}

//----------------------------------------------------------------------------
// Get primitive variable from a conserved variable
//----------------------------------------------------------------------------

void Conservation_Law::conserved_to_primitive(const Vector<double>& U, Vector<double>& W) const {
/*
  Assert(W.size() == 5,
         ExcDimensionMismatch(W.size(), 5) );
  Assert(U.size() == 5,
         ExcDimensionMismatch(U.size(), 5) );
*/
    W[0] = U[0];
    W[1] = U[1]/U[0];
    W[2] = U[2]/U[0];
    W[3] = (gamma()-1.0)*(U[3] - 0.5*((U[1]*U[1] + U[2]*U[2])/U[0]));
    W[4] = U[4];
    //W[3] = p/(W[0]*R());
}

void Conservation_Law::conserved_to_primitive(const Vector<double>& U, double gamma, double p_inf, Vector<double>& W) const {
/*
  Assert(W.size() == 5,
         ExcDimensionMismatch(W.size(), 5) );
  Assert(U.size() == 5,
         ExcDimensionMismatch(U.size(), 5) );
*/
    W[0] = U[0];
    W[1] = U[1]/U[0];
    W[2] = U[2]/U[0];
    W[3] = (gamma -1.0)*(U[3] - 0.5*((U[1]*U[1] + U[2]*U[2])/U[0]))  - gamma*p_inf;
    W[4] = U[4];
    //W[3] = p/(W[0]*R());
}

//----------------------------------------------------------------------------
// Compute purely convective part of the flux
//----------------------------------------------------------------------------

void Conservation_Law::compute_pure_convective_flux(const Vector<double>& U, const double nx, const double ny, double x, double y, Vector<double>& Flux) {

	double p = (gamma()-1.0)*(U[3] - 0.5*(U[1]*U[1] + U[2]*U[2])/U[0]);


	ErrNegativePressureDensity(U[0], p, x, y);

	double u = U[1]/U[0];
	double v = U[2]/U[0];

	double F_c[5]; double G_c[5];

	F_c[0] = U[1];                 // rho*u
	F_c[1] = U[0]*u*u + p ;        // rho*u^2 + p
	F_c[2] = U[0]*u*v;             // rho*u*v
	F_c[3] = u*(p + U[3]);         // u*(E + p)
	F_c[4] = u*U[4];                 // u*phi

	G_c[0] = U[2];                 // rho*v
	G_c[1] = F_c[2] ;              // rho*u*v
	G_c[2] = U[0]*v*v + p ;        // rho*v^2 + p
	G_c[3] = v*(p + U[3]);         // v*(E + p)
	G_c[4] = v*U[4];                 // rho*v

	for (unsigned int n_components = 0; n_components < 5; ++n_components) {
		Flux[n_components] = nx*F_c[n_components] + ny*G_c[n_components];
	}

}


//----------------------------------------------------------------------------
// Compute purely viscous part of the flux
//----------------------------------------------------------------------------

void Conservation_Law::compute_pure_viscous_flux(
		const Vector<double>& U,
		const FullMatrix<double>& grad_U,
		const double gamma_L, const double p_inf_L,
		const double nx, const double ny,
		Vector<double>& Flux) {

	static const double r2_3 = 2./3, r4_3 = 4./3;

	double u = U[1]/U[0];
	double v = U[2]/U[0];

	double mu = viscosity(U);
	double k = heat_conductivity(U);


	// Viscous flux (for derivations see the sympy notebook)

	double tau_xx = mu*(r4_3*grad_U(1,0) - r2_3*grad_U(2,1)) ;

	double tau_xy = mu*(grad_U(1,1) + grad_U(2,0));

	double tau_yy = mu*(r4_3*grad_U(2,1) - r2_3*grad_U(1,0)) ;

	double q_x = -k*grad_U(3,0);

	double q_y = -k*grad_U(3,1);

	double F_v[4]; double G_v[4];

	F_v[0] =   0.0;
	F_v[1] = - tau_xx;
	F_v[2] = - tau_xy;
	F_v[3] = - (u*tau_xx + v*tau_xy - q_x);


	G_v[0] =   0.0;
	G_v[1] = - tau_xy;
	G_v[2] = - tau_yy;
	G_v[3] = - (u*tau_xy + v*tau_yy - q_y);

	for (unsigned int n_components = 0; n_components < 4; ++n_components) {
		Flux[n_components] = nx*F_v[n_components] + ny*G_v[n_components];
	}
}

//----------------------------------------------------------------------------
// Find the temperature gradient at a point
//----------------------------------------------------------------------------

void Conservation_Law::temperature_gradient(
		const Vector<double>& U,
		const FullMatrix<double>& grad_U,
		double* grad_T) {

	grad_T[0] = (gamma() - 1.0)*(-2*(U(1)*grad_U(1,0) + U(2)*grad_U(2,0))*U(0) +
			(U(1)*U(1) + U(2)*U(2))*grad_U(0,0) +
			(-2*U(3)*U(0) + U(1)*U(1) + U(2)*U(2))*grad_U(0,0) +
			2*U(0)*U(0)*grad_U(3,0))/(2*U(0)*U(0)*U(0));

	grad_T[1] = (gamma() - 1.0)*(-2*(U(1)*grad_U(1,1) + U(2)*grad_U(2,1))*U(0) +
			(U(1)*U(1) + U(2)*U(2))*grad_U(0,1) +
			(-2*U(3)*U(0) + U(1)*U(1) + U(2)*U(2))*grad_U(0,1) +
			2*U(0)*U(0)*grad_U(3,1))/(2*U(0)*U(0)*U(0));
}

//----------------------------------------------------------------------------
// Compute the fluxes F(U, grad_U) and G(U, grad_U)
//----------------------------------------------------------------------------

void Conservation_Law::compute_x_y_flux(
		const Vector<double>& U,
		const FullMatrix<double>& grad_U,
		double x, double y,
		Vector<double>& F,
		Vector<double>& G) {

	static const double r2_3 = 2./3, r4_3 = 4./3;

	double p = (gamma()-1.0)*(U[3] - 0.5*(U[1]*U[1] + U[2]*U[2])/U[0]);

	ErrNegativePressureDensity(U[0], p, x, y);

	// Convective flux

	double u = U[1]/U[0];
	double v = U[2]/U[0];

	double mu = viscosity(U);
	double k = heat_conductivity(U);

	F[0] = U[1];                 // rho*u
	F[1] = U[0]*u*u + p ;        // rho*u^2 + p
	F[2] = U[0]*u*v;             // rho*u*v
	F[3] = u*(p + U[3]);         // u*(E + p)


	G[0] = U[2];                 // rho*v
	G[1] = F[2] ;                // rho*u*v
	G[2] = U[0]*v*v + p ;        // rho*v^2 + p
	G[3] = v*(p + U[3]);         // v*(E + p)

	// Viscous flux (for derivations see the sympy notebook)

	double tau_xx = mu*(r4_3*grad_U(1,0) - r2_3*grad_U(2,1)) ;

	double tau_xy = mu*(grad_U(1,1) + grad_U(2,0));

	double tau_yy = mu*(r4_3*grad_U(2,1) - r2_3*grad_U(1,0)) ;

	double q_x = -k*grad_U(3,0);

	double q_y = -k*grad_U(3,1);

/*
	// comments start
	double rho = U[0];
	double M_x = U[1];
	double M_y = U[2];
	double E   = U[3];
	double rho_x = grad_U[0][0]; double rho_y = grad_U[0][1];
	double M_x_x = grad_U[1][0]; double M_x_y = grad_U[1][1];
	double M_y_x = grad_U[2][0]; double M_y_y = grad_U[2][1];
	double   E_x = grad_U[3][0]; double   E_y = grad_U[3][1];

	double tau_xx = 2.0*mu*(-2.0*M_x*rho_x + M_y*rho_y + 2*rho*M_x_x - rho*M_y_y)/(3.0*rho*rho);
	double tau_xy = -mu*(-(M_x_y + M_y_x)*rho + M_x*rho_y + M_y*rho_x)/(rho*rho);
	double tau_yy = 2*mu*(M_x*rho_x - 2*M_y*rho_y - rho*M_x_x + 2*rho*M_y_y)/(3.0*rho*rho);


	double q_x = -heat_conductivity(T)*(gamma() - 1)*(-2*(M_x*M_x_x + M_y*M_y_x)*rho + (M_x*M_x + M_y*M_y)*rho_x + (-2.0*E*rho + M_x*M_x + M_y*M_y)*rho_x + 2.0*rho*rho*E_x)/(2.0*R()*rho*rho*rho);
	double q_y = -heat_conductivity(T)*(gamma() - 1)*(-2*(M_x*M_x_y + M_y*M_y_y)*rho + (M_x*M_x + M_y*M_y)*rho_y + (-2.0*E*rho + M_x*M_x + M_y*M_y)*rho_y + 2.0*rho*rho*E_y)/(2.0*R()*rho*rho*rho);
	// comments end
*/
	F[1] = F[1] - tau_xx;
	F[2] = F[2] - tau_xy;
	F[3] = F[3] - (u*tau_xx + v*tau_xy - q_x);

	G[1] = G[1] - tau_xy;
	G[2] = G[2] - tau_yy;
	G[3] = G[3] - (u*tau_xy + v*tau_yy - q_y);
}

//----------------------------------------------------------------------------
// Compute the normal component of the flux:
// F(U, grad_U)*nx + G(U, grad_U)*ny
//----------------------------------------------------------------------------

void Conservation_Law::compute_physical_flux(
		const Vector<double>& U,
		const FullMatrix<double>& grad_U,
		double nx, double ny,
		double x, double y,
		Vector<double>& Flux) {

	Vector<double> F(5); Vector<double> G(5);

	compute_x_y_flux(U, grad_U, x, y, F, G);

	for (unsigned int i = 0; i < 4; ++i) {
		Flux[i] = nx*F[i] + ny*G[i];
	}
}


//----------------------------------------------------------------------------
// Local Lax–Friedrichs Riemann solver for convective part of NS equations
//----------------------------------------------------------------------------

void Conservation_Law::LLF_riemann_solver(
		const Vector<double>& UL,
		const Vector<double>& UR,
		double nx, double ny,
		double x, double y, Vector<double>& Flux) {

	Vector<double> FL(5); Vector<double> FR(5);

	compute_pure_convective_flux(UL, nx, ny, x, y, FL);
	compute_pure_convective_flux(UR, nx, ny, x, y, FR);

	// Maximum convective speed

	double a_l = speed_of_sound(UL);
	double a_r = speed_of_sound(UR);

	double s_l = std::abs(UL[1]*nx + UL[2]*ny)/UL[0] + a_l;
	double s_r = std::abs(UR[1]*nx + UR[2]*ny)/UR[0] + a_r;

	double s_max = std::max(s_l, s_r);

	for (unsigned int i = 0; i < 4; ++i) {
		Flux[i] = 0.5*(FL[i] + FR[i] - s_max*(UR[i] - UL[i]));
	}
}

//----------------------------------------------------------------------------
// Apply transformation and inverse transformation in normal direction
//----------------------------------------------------------------------------

void multiply_with_rotation_matrix(const Vector<double>& U, double nx, double ny, Vector<double>& U_hat) {

	U_hat[0] = U[0];
    U_hat[1] = U[1]*nx + U[2]*ny;
    U_hat[2] = -U[1]*ny + U[2]*nx;
    U_hat[3] = U[3];
	U_hat[4] = U[4];
}

void multiply_with_inv_rotation_matrix(const Vector<double>& U, double nx, double ny, Vector<double>& U_hat) {

	U_hat[0] = U[0];
    U_hat[1] = U[1]*nx - U[2]*ny;
    U_hat[2] = U[1]*ny + U[2]*nx;
    U_hat[3] = U[3];
	U_hat[4] = U[4];
	U_hat[5] = U[5];
}

//----------------------------------------------------------------------------
// HLLC Riemann solver for convective part of NS equations
//----------------------------------------------------------------------------
/*
void Conservation_Law::HLLC_riemann_solver(
		const Vector<double>& UL,
		const Vector<double>& UR,
		double nx, double ny, double x, double y, Vector<double>& Flux) {

	Vector<double> U_hatL(5); Vector<double> U_hatR(5); Vector<double> F_hat(5);

    multiply_with_rotation_matrix(UL, nx, ny, U_hatL);
    multiply_with_rotation_matrix(UR, nx, ny, U_hatR);

    double g3 = (gamma() + 1.0)/(2*gamma());

	// Extract states
	double rho_L = UL[0]; double u_L = U_hatL[1]/UL[0]; double v_L = U_hatL[2]/UL[0]; double p_L = Pressure(UL);
	double rho_R = UR[0]; double u_R = U_hatR[1]/UR[0]; double v_R = U_hatR[2]/UR[0]; double p_R = Pressure(UR);
	double phi_L = UL[4];	double phi_R = UR[4];

	ErrNegativePressureDensity(rho_L, p_L, x, y);
	ErrNegativePressureDensity(rho_R, p_R, x, y);

	double a_L = std::sqrt(gamma()*p_L/rho_L);
	double a_R = std::sqrt(gamma()*p_R/rho_R);

	// Pressure Estimates

	double rho_bar  = 0.5*(rho_L + rho_R); double a_bar = 0.5*(a_L + a_R);
	double p_pvrs = 0.5*(p_L + p_R) - 0.5*(u_R - u_L)*rho_bar*a_bar;
	double p_star = std::max(0.0, p_pvrs);

	// Wave speed estimates

	double q_L, q_R;

	if (p_star <= p_L) {
		q_L = 1.0;
	}
	else {
		q_L = std::sqrt(1.0 + g3*((p_star/p_L) - 1.0));
	}

	if (p_star <= p_R) {
		q_R = 1.0;
	}
	else {
		q_R = std::sqrt(1.0 + g3*((p_star/p_R) - 1.0));
	}

	double S_L = u_L - a_L*q_L; double S_R = u_R + a_R*q_R;

	double S_star = (p_R - p_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R))
			/(rho_L*(S_L - u_L) - rho_R*(S_R- u_R));

	// HLLC Flux

	double FL[6]; double FR[6];

	FL[0] = U_hatL[1];
	FL[1] = UL[0]*u_L*u_L + p_L;
	FL[2] = UL[0]*u_L*v_L;
	FL[3] = u_L*(p_L + UL[3]);
	FL[4] = u_L*U_hatL[4];

	FR[0] = U_hatR[1];
	FR[1] = rho_R*u_R*u_R + p_R;
	FR[2] = rho_R*u_R*v_R;
	FR[3] = u_R*(p_R + UR[3]);
	FR[4] = u_R*U_hatR[4];

    if (0.0 <= S_L) {

        F_hat[0] = FL[0];
        F_hat[1] = FL[1];
        F_hat[2] = FL[2];
        F_hat[3] = FL[3];
        F_hat[4] = FL[4];
		F_hat[5] = 0.0;
    }

    else if (S_L <= 0.0 && 0.0 <= S_star) {

        double U_starL[5];

        U_starL[0] = rho_L*((S_L - u_L)/(S_L - S_star));
        U_starL[1] = rho_L*((S_L - u_L)/(S_L - S_star))*S_star;
        U_starL[2] = rho_L*((S_L - u_L)/(S_L - S_star))*v_L;
        U_starL[3] = rho_L*((S_L - u_L)/(S_L - S_star))*(U_hatL[3]/rho_L + (S_star - u_L)*(S_star + p_L/(rho_L*(S_L - u_L))));
        U_starL[4] = phi_L*((S_L - u_L)/(S_L - S_star));

        F_hat[0] = FL[0] + S_L*(U_starL[0] - U_hatL[0]);
        F_hat[1] = FL[1] + S_L*(U_starL[1] - U_hatL[1]);
        F_hat[2] = FL[2] + S_L*(U_starL[2] - U_hatL[2]);
        F_hat[3] = FL[3] + S_L*(U_starL[3] - U_hatL[3]);
        F_hat[4] = FL[4] + S_L*(U_starL[4] - U_hatL[4]);
		F_hat[5] = u_L + S_L*( ( (S_L - u_L)/(S_L - S_star) ) - 1.0 ) ;
    }


    else if (S_star <= 0.0 && 0.0 <= S_R) {

    	double U_starR[5];

        U_starR[0] = rho_R*((S_R - u_R)/(S_R - S_star));
        U_starR[1] = rho_R*((S_R - u_R)/(S_R - S_star))*S_star;
        U_starR[2] = rho_R*((S_R - u_R)/(S_R - S_star))*v_R;
        U_starR[3] = rho_R*((S_R - u_R)/(S_R - S_star))*(U_hatR[3]/rho_R + (S_star - u_R)*(S_star + p_R/(rho_R*(S_R - u_R))));
        U_starR[4] = phi_R*((S_R - u_R)/(S_R - S_star));

        F_hat[0] = FR[0] + S_R*(U_starR[0] - U_hatR[0]);
        F_hat[1] = FR[1] + S_R*(U_starR[1] - U_hatR[1]);
        F_hat[2] = FR[2] + S_R*(U_starR[2] - U_hatR[2]);
        F_hat[3] = FR[3] + S_R*(U_starR[3] - U_hatR[3]);
        F_hat[4] = FR[4] + S_R*(U_starR[4] - U_hatR[4]);
		F_hat[5] = u_R + S_R*( ( (S_R - u_R)/(S_R - S_star) ) - 1.0 ) ;
    }

    else if (0.0 >= S_R) {

        F_hat[0] = FR[0];
        F_hat[1] = FR[1];
        F_hat[2] = FR[2];
        F_hat[3] = FR[3];
        F_hat[4] = FR[4];
		F_hat[5] = 0.0;
    }

    multiply_with_inv_rotation_matrix(F_hat, nx, ny, Flux);
}
*/

void Conservation_Law::FluxFunction(const Vector<double>& W, double gamma_l, double P_Infty_l, double nx, double ny, Vector<double>& FV) {
	double un; // normal velocity
	un = (W[1]*nx + W[2]*ny) ;
        FV[0] = W[0]*un ;
        FV[1] = W[0]*W[1]*un + W[3]*nx ;
        FV[2] = W[0]*W[2]*un + W[3]*ny ;
        FV[3] = ( gamma_l*(W[3] + P_Infty_l)/(gamma_l - 1.0) + 0.5*W[0]*(W[1]*W[1] + W[2]*W[2]) )*un ;
      FV[4] = W[4]*un ;
}

void Conservation_Law::HLLC_riemann_solver(
		const Vector<double>& WL,
		const Vector<double>& WR,
		double gamma_L, double gamma_R,
		double p_inf_L, double p_inf_R,
		double nx, double ny, double x, double y, Vector<double>& Flux) {   

	unsigned int size_ = 5, size_F = 6;
    
    int i ;
	double rho_L, rho_R, u_L, u_R, v_L, v_R, P_L, P_R, c_L, c_R, E_L, E_R ;
	double un_L, un_R, ut_L, ut_R ;
	double un, ut ;
	double S_L, S_R, S_star;
	Vector<double> UL(size_), UR(size_), UL_star(size_), UR_star(size_), FL(size_F), FR(size_F), FL_star(size_F), FR_star(size_F);

//	conserved_to_primitive(UL, gamma_L, p_inf_L, WL);
//	conserved_to_primitive(UR, gamma_R, p_inf_R, WR);
	
	rho_L = WL[0] ; u_L = WL[1] ; v_L = WL[2] ; P_L = WL[3] ;
    rho_R = WR[0] ; u_R = WR[1] ; v_R = WR[2] ; P_R = WR[3] ;

//	ErrNegativePressureDensity(rho_L, P_L+p_inf_L, x, y);
//	ErrNegativePressureDensity(rho_R, P_R+p_inf_R, x, y);

	un_L = u_L*nx + v_L*ny ; ut_L = -u_L*ny + v_L*nx ;
	un_R = u_R*nx + v_R*ny ; ut_R = -u_R*ny + v_R*nx ;

	if( (P_L < -p_inf_L) || (P_R < -p_inf_R) ) {
		std::cout << "\n Imaginary speed of sound " << std::endl ;
		std::cout << " Density : " << rho_L <<  "\t" << rho_R << std::endl ;
		std::cout << " Normal Velocity: " << un_L <<  "\t" << un_R << std::endl ;
		std::cout << " Pressure: " << P_L <<  "\t" << P_R << std::endl ;
		std::cout << " Phi: " << WL[4] <<  "\t" << WR[4] << std::endl ;
		std::cout << " Gamma:    " << gamma_L <<  "\t" << gamma_R << std::endl ;
		std::cout << " P_Infty:  " << p_inf_L <<  "\t" << p_inf_R << std::endl ;
		std::cout << " x:  " << x <<  "\ty" << y << std::endl ;
		std::cout << " nx:  " << nx <<  "\tny" << ny << std::endl ;		
		std::exit(1);
	}

    c_L = sqrt(gamma_L*(P_L + p_inf_L)/rho_L) ; c_R = sqrt(gamma_R*(P_R + p_inf_R)/rho_R) ;

	E_L = (P_L + gamma_L*p_inf_L)/(gamma_L - 1.0) + 0.5*rho_L*u_L*u_L + 0.5*rho_L*v_L*v_L ;
	E_R = (P_R + gamma_R*p_inf_R)/(gamma_R - 1.0) + 0.5*rho_R*u_R*u_R + 0.5*rho_R*v_R*v_R ;

	S_L = std::min((un_R - c_R), (un_L - c_L)) ; 
	S_R = std::max((un_L + c_L), (un_R + c_R)) ;
	S_star = ( P_R - P_L + rho_L*un_L*(S_L-un_L) - rho_R*un_R*(S_R - un_R) )/(rho_L*(S_L - un_L) - rho_R*(S_R - un_R)) ;
	//P_Star = P_L + rho_L*(S_star - un_L)*(S_L - un_L) ;

	// Now compute the left right and starred fluxes for HLLC.
	FluxFunction(WL,gamma_L,p_inf_L,nx,ny, FL) ; FluxFunction(WR,gamma_R,p_inf_R,nx,ny, FR);
	FL[5] = un_L ; FR[5] = un_R ;

	UL_star[0] = rho_L*(S_L - un_L)/(S_L - S_star) ;  
	un = S_star ; ut = ut_L ;
	UL_star[1] = UL_star[0]*(un*nx - ut*ny) ;  
	UL_star[2] = UL_star[0]*(un*ny + ut*nx) ;
	UL_star[3] = UL_star[0]*( (E_L/rho_L) + (S_star - un_L)*(S_star + P_L/(rho_L*(S_L - un_L)) ) )	;  
    UL_star[4] = WL[4]*(S_L - un_L)/(S_L - S_star);	

	UR_star[0] = rho_R*(S_R - un_R)/(S_R - S_star) ;
	un = S_star ; ut = ut_R ;
	UR_star[1] = UR_star[0]*(un*nx - ut*ny) ;
	UR_star[2] = UR_star[0]*(un*ny + ut*nx) ;
	UR_star[3] = UR_star[0]*( (E_R/rho_R) + (S_star - un_R)*(S_star + P_R/(rho_R*(S_R - un_R)) ) ) ;
    UR_star[4] = WR[4]*(S_R - un_R)/(S_R - S_star) ;

	primitive_to_conserved(WL,gamma_L,p_inf_L,UL) ;
	for(i = 0 ; i < 5 ; i++) FL_star[i] = FL[i] + S_L*(UL_star[i] - UL[i]) ; 
    FL_star[5] = un_L + S_L*( ((S_L - un_L)/(S_L - S_star)) - 1.0 ) ;
	primitive_to_conserved(WR,gamma_R,p_inf_R,UR) ;
	for(i = 0 ; i < 5 ; i++) FR_star[i] = FR[i] + S_R*(UR_star[i] - UR[i]) ; 
    FR_star[5] = un_R + S_R*( ((S_R - un_R)/(S_R - S_star)) - 1.0 ) ;

	if( S_L > 0.0 ) {
		for(i = 0 ; i < 6 ; i++) Flux[i] = FL[i] ; 
	} else if((S_star >= 0.0) && (S_L < 0.0)) {
		for(i = 0 ; i < 6 ; i++) Flux[i] = FL_star[i] ; 
	} else if((S_star < 0.0) && (S_R >= 0.0)) {
		for(i = 0 ; i < 6 ; i++) Flux[i] = FR_star[i] ; 
	} else if(S_R < 0.0) {
		for(i = 0 ; i < 6 ; i++) Flux[i] = FR[i] ; 
	}

}


//----------------------------------------------------------------------------
// Rotated HLLC Riemann solver for convective part of NS equations
//----------------------------------------------------------------------------


//Vector<double> rotated_HLLC_riemann_solver(Vector<double> UL, Vector<double> UR, double nx, double ny, Point<2> P, bool boundary) {
void Conservation_Law::rotated_HLLC_riemann_solver(
		const Vector<double>& WL,
		const Vector<double>& WR,
		double gamma_L, double gamma_R,
		double p_inf_L, double p_inf_R,
		double nx, double ny, double x, double y, Vector<double>& Flux) { 

    Vector<double> flux1(6);
    Vector<double> flux2(6);

    int i ;
	double alpha1, alpha2, n1x, n1y, n2x, n2y, u_L, u_R, v_L, v_R, du, dv, dq ;

	u_L = WL[1] ; v_L = WL[2] ; 
    u_R = WR[1] ; v_R = WR[2] ; 
	du = u_R - u_L ; dv = v_R - v_L ; 
	dq = sqrt(du*du + dv*dv) ;
	
	if(dq < 1.0E-10) { n1x = nx ; n1y = ny ; }
	else { n1x = du/dq ; n1y = dv/dq ; }

	alpha1 = (n1x*nx + n1y*ny) ;
	if(alpha1 < 0) { n1x = -n1x ; n1y = -n1y ; alpha1 = -alpha1 ; }
	n2x = -n1y ; n2y = n1x ;
	alpha2 = (n2x*nx + n2y*ny) ; 
	if(alpha2 < 0) { n2x = -n2x ; n2y = -n2y ; alpha2 = -alpha2 ; }

	HLLC_riemann_solver(WL, WR, gamma_L, gamma_R, p_inf_L, p_inf_R, n1x, n1y, x, y, flux1);
	HLLC_riemann_solver(WL, WR, gamma_L, gamma_R, p_inf_L, p_inf_R, n2x, n2y, x, y, flux2);

	for(i = 0 ; i < 6 ; i++) Flux[i] = alpha1*flux1[i] + alpha2*flux2[i] ;
//	for(i = 0 ; i < 4 ; i++) Flux[i] = alpha1*flux1[i] + alpha2*flux2[i] ;
}


//----------------------------------------------------------------------------
// Riemann solver for viscous part of NS equations
//----------------------------------------------------------------------------

void Conservation_Law::viscous_riemann_solver(
		const Vector<double>& UL,
		const Vector<double>& UR,
		const FullMatrix<double>& grad_UL,
		const FullMatrix<double>& grad_UR,
		double gamma_L, double gamma_R,
		double p_inf_L, double p_inf_R,
		double h, double nx, double ny, unsigned int order, Vector<double>& Flux) {


	Vector<double> FL(4); Vector<double> FR(4);

	compute_pure_viscous_flux(UL, grad_UL, gamma_L, p_inf_L, nx, ny, FL);
	compute_pure_viscous_flux(UR, grad_UR, gamma_R, p_inf_R, nx, ny, FR);

	for (unsigned int i = 0; i < 4; ++i) {
		Flux[i] = 0.5*(FL[i] + FR[i]);
	}

	Flux[4] = 0.0; Flux[5] = 0.0;
}

//----------------------------------------------------------------------------
// Combine the numerical viscous and convective fluxes to get the final flux
//----------------------------------------------------------------------------

Vector<double> Conservation_Law::numerical_flux(
		const Vector<double>& WL,
		const Vector<double>& WR,
		const FullMatrix<double>& grad_WL,
		const FullMatrix<double>& grad_WR,
		double gamma_L, double gamma_R,
		double p_inf_L, double p_inf_R,
		double h, double nx, double ny, double x, double y, unsigned int order) {


	Vector<double> Flux(6);

	Vector<double> F_c(6); Vector<double> F_v(6);

//	LLF_riemann_solver(UL, UR, nx, ny, x, y, F_c);
	rotated_HLLC_riemann_solver(WL, WR, gamma_L, gamma_R, p_inf_L, p_inf_R, nx, ny, x, y, F_c);
/*
	viscous_riemann_solver(WL, WR, grad_WL, grad_WR, gamma_L, gamma_R, p_inf_L, p_inf_R, h, nx, ny, order, F_v);

	for (unsigned int n_component = 0; n_component < 6; ++n_component) {
		Flux[n_component] = F_c[n_component] + F_v[n_component];
	}

	return Flux;
*/
	return F_c;
}

void Conservation_Law::axisymmetric_source(Point<2> center, const Vector<double>& W, Vector<double>& source) {
    double y = center(1);
    double phi = W(4);
    double gamma_Fluid = 1.0 + (GAMMA_L-1.0)*(GAMMA_G-1.0)/((1.0-phi)*(GAMMA_L -1.0) + phi*(GAMMA_G-1.0));
    double P_inf = ((gamma_Fluid -1.0)/gamma_Fluid)*( GAMMA_L*P_INF_L*phi/(GAMMA_L - 1.0) + GAMMA_G*P_INF_G*(1.0 - phi)/(GAMMA_G - 1.0) );
    double rho = W(0);
    double u = W(1);
    double v = W(2);
    double p = W(3);
    double E = (p + gamma_Fluid*P_inf)/(gamma_Fluid - 1.0) + 0.5*rho*(u*u + v*v);
    double miy = -1.0/y;

    source(0) = miy*rho*v;
    source(1) = miy*rho*u*v;
    source(2) = miy*rho*v*v;
    source(3) = miy*(E+p)*v;
    source(4) = 0.0;
}

