/*
 * pde.c
 *      Author: sunder
 */

#include "hype.h"

//----------------------------------------------------------------------------
// Convert a conserved variable to primitive variable 
//----------------------------------------------------------------------------

void PDECons2Prim(const PetscReal *Q, PetscReal *V) {
    
    PetscReal phi = Q[4]; 
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-phi)*(g1 -1.0) + phi*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*phi/(g1 - 1.0) + g2*p2*(1.0 - phi)/(g2 - 1.0) );
    PetscReal irho = 1.0/Q[0]; 
    
    V[0] = Q[0];
    V[1] = irho*Q[1];
    V[2] = irho*Q[2];
    V[3] = (g -1.0)*( Q[3] - 0.5*irho*(Q[1]*Q[1] + Q[2]*Q[2]) )  - g*P_inf;
    V[4] = phi;
    
}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable 
//----------------------------------------------------------------------------

void PDEPrim2Cons(const PetscReal *V, PetscReal *Q) {

    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-V[4])*(g1 -1.0) + V[4]*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*V[4]/(g1 - 1.0) + g2*p2*(1.0 - V[4])/(g2 - 1.0) );

    PetscReal e = (V[3] + g*P_inf)/(g - 1.0);
    PetscReal k = 0.5*V[0]*(V[1]*V[1] + V[2]*V[2]);

    Q[0] = V[0];
    Q[1] = V[0]*V[1];
    Q[2] = V[0]*V[2];
    Q[3] = k + e;
    Q[4] = V[4];
}

//----------------------------------------------------------------------------
// Conservative flux components F in the given normal direction 
//----------------------------------------------------------------------------

PetscReal PDEFluxPrim(const PetscReal *V, 
                            const PetscReal nx, const PetscReal ny,
                            const PetscReal x,  const PetscReal y,
                            PetscReal *F) {

    PetscReal phi = V[4]; 
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-phi)*(g1 -1.0) + phi*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*phi/(g1 - 1.0) + g2*p2*(1.0 - phi)/(g2 - 1.0) );
    PetscReal rho = V[0]; 
    PetscReal u = V[1];
    PetscReal v = V[2];
    PetscReal p = V[3];
    PetscReal E = (p + g*P_inf)/(g - 1.0) + 0.5*rho*(u*u + v*v); 
    PetscReal un = u*nx + v*ny; 
    PetscReal rhoun = rho*un; 
    
    // Check if the input state is physically admissible 

    if (rho < rho_floor) {
        printf("Negative density = %f\n", rho);
        printf("At x = %f, y = %f\n", x, y);  
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    if ((p + P_inf)  < prs_floor) {
        printf("Negative pressure, p + p_inf = %f\n", p + P_inf);
        printf("At x = %f, y = %f\n", x, y);  
        MPI_Abort(PETSC_COMM_WORLD, 1);
    }

    // Now find the fluxes 

    F[0] = rhoun;
    F[1] = rhoun*u + p*nx;
    F[2] = rhoun*v + p*ny;
    F[3] = un*(E + p);
    F[4] = phi*un;

    // Also obtain the maximum eigen value 

    PetscReal s_max = PetscAbsReal(un) + PetscSqrtReal(g*(p + P_inf)/rho);

    return s_max;
}

//----------------------------------------------------------------------------
// Viscous part of the flux in the normal direction 
//----------------------------------------------------------------------------

PetscReal PDEViscFluxPrim(PetscReal y, const PetscReal* V, const PetscReal grad_V[nVar][DIM], PetscReal nx, PetscReal ny, PetscReal* F) {
    
    PetscReal r2_3 = 2./3.; PetscReal r4_3 = 4./3.;

    // Find the phase fractions

    PetscReal phi = V[4];

    // Effective viscosity

    PetscReal mu = phi*mu1 + (1.0-phi)*mu2;

    PetscReal u = V[1];
    PetscReal v = V[2];
    PetscReal u_x = grad_V[1][0];
    PetscReal u_y = grad_V[1][1];
    PetscReal v_x = grad_V[2][0];
    PetscReal v_y = grad_V[2][1];
    PetscReal v_over_y;

    if (PetscAbsReal(y) < small_num)
        v_over_y = v_y; 
    else
        v_over_y = v/y; 

    PetscReal div_v = u_x + v_y + v_over_y;

    // Stress tensor

    PetscReal tau_xx = mu*(2.0*u_x - r2_3*div_v);
    PetscReal tau_xy = mu*(u_y + v_x);
    PetscReal tau_yy = mu*(2.0*v_y - r2_3*div_v);

    F[0] = 0.0;
    F[1] = -nx*tau_xx - ny*tau_xy;
    F[2] = -nx*tau_xy - ny*tau_yy;
    F[3] = -nx*(u*tau_xx + v*tau_xy) - ny*(u*tau_xy + v*tau_yy);
    F[4] = 0.0;
    F[5] = 0.0;

    return r4_3*mu/V[0];
}

//----------------------------------------------------------------------------
// Source terms
//----------------------------------------------------------------------------

void PDESource(PetscReal y, const PetscReal *V, const PetscReal *grad_V_x, const PetscReal *grad_V_y, PetscReal *S) {

    PetscReal miy = -1.0/y;
    PetscReal iy = 1.0/y;

    PetscReal phi = V[4];
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-phi)*(g1 -1.0) + phi*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*phi/(g1 - 1.0) + g2*p2*(1.0 - phi)/(g2 - 1.0) );
    PetscReal rho = V[0];
    PetscReal u = V[1];
    PetscReal v = V[2];
    PetscReal p = V[3];
    PetscReal E = (p + g*P_inf)/(g - 1.0) + 0.5*rho*(u*u + v*v);

    PetscReal r2_3 = 2./3.;
    PetscReal mu = phi*mu1 + (1.0-phi)*mu2;

    // Viscous Source Terms

    PetscReal u_x = grad_V_x[1];
    PetscReal u_y = grad_V_y[1];
    PetscReal v_x = grad_V_x[2];
    PetscReal v_y = grad_V_y[2];
    PetscReal div_v = u_x + v_y + v*iy;

    // Stress tensor

    PetscReal tau_zz = mu*(2.0*v*iy - r2_3*div_v);
    PetscReal tau_xy = mu*(u_y + v_x) ;
    PetscReal tau_yy = mu*(2.0*v_y - r2_3*div_v);

    S[0] = miy*rho*v;
    S[1] = miy*rho*u*v + iy*tau_xy;
    S[2] = miy*rho*v*v + iy*(tau_yy - tau_zz);
    S[3] = miy*v*(E+p) + iy*(v*tau_yy + u*tau_xy);
    S[4] = 0.0;

}

//----------------------------------------------------------------------------
// Check physical admissibility of the input state 
//----------------------------------------------------------------------------

PetscBool PDECheckPADPrim(const PetscReal *V) {
    
    PetscBool PAD = PETSC_TRUE;
    
    PetscReal phi = V[4]; 
    PetscReal g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-phi)*(g1 -1.0) + phi*(g2-1.0));
    PetscReal P_inf = ((g -1.0)/g)*( g1*p1*phi/(g1 - 1.0) + g2*p2*(1.0 - phi)/(g2 - 1.0) );
    PetscReal rho = V[0]; 
    PetscReal p = V[3];

    // Check if the input state is physically admissible 


    if (rho < rho_floor) {    
        PAD = PETSC_FALSE;
    }

    if ((p + P_inf)  < prs_floor) {
        PAD = PETSC_FALSE;
    }
    
    if (phi < 0.0 || phi > 1.0) {
        PAD = PETSC_FALSE;
    }
    
    return PAD; 
}

//----------------------------------------------------------------------------
// HLLC Riemann solver 
//----------------------------------------------------------------------------

PetscReal HLLCRiemannSolver(const PetscReal *VL, const PetscReal *VR,
               const PetscReal nx, const PetscReal ny,
               const PetscReal x,  const PetscReal y,
               PetscReal *Flux) {
    
    
    PetscReal QL[nVar], QR[nVar]; 
    PetscReal FL[6], FR[6]; 
    const PetscInt size_ = 5, size_F = 6;
    PetscInt i;
    PetscReal rho_L, rho_R, u_L, u_R, v_L, v_R, P_L, P_R, c_L, c_R, E_L, E_R ;
    PetscReal un_L, un_R, ut_L, ut_R ;
    PetscReal un, ut ;
    PetscReal S_L, S_R, S_star;
    PetscReal UL_star[size_], UR_star[size_], FL_star[size_F], FR_star[size_F];
    
    PDEPrim2Cons(VL, QL); PDEPrim2Cons(VR, QR);
    
    PetscReal s_max_l = PDEFluxPrim(VL, nx, ny, x, y, FL); 
    PetscReal s_max_r = PDEFluxPrim(VR, nx, ny, x, y, FR);
    
    PetscReal gamma_L = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-VL[4])*(g1 -1.0) + VL[4]*(g2-1.0));
	PetscReal p_inf_L = ((gamma_L -1.0)/gamma_L)*( g1*p1*VL[4]/(g1 - 1.0) + g2*p2*(1.0 - VL[4])/(g2 - 1.0) );
    PetscReal gamma_R = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-VR[4])*(g1 -1.0) + VR[4]*(g2-1.0));
    PetscReal p_inf_R = ((gamma_R -1.0)/gamma_R)*( g1*p1*VR[4]/(g1 - 1.0) + g2*p2*(1.0 - VR[4])/(g2 - 1.0) );
    
    PetscReal s_max = PetscMax(s_max_l, s_max_r);
    
    rho_L = VL[0] ; u_L = VL[1] ; v_L = VL[2] ; P_L = VL[3];
    rho_R = VR[0] ; u_R = VR[1] ; v_R = VR[2] ; P_R = VR[3];

    un_L = u_L*nx + v_L*ny ; ut_L = -u_L*ny + v_L*nx;
    un_R = u_R*nx + v_R*ny ; ut_R = -u_R*ny + v_R*nx;

    c_L = PetscSqrtReal(gamma_L*(P_L + p_inf_L)/rho_L) ; c_R = PetscSqrtReal(gamma_R*(P_R + p_inf_R)/rho_R) ;

    E_L = (P_L + gamma_L*p_inf_L)/(gamma_L - 1.0) + 0.5*rho_L*u_L*u_L + 0.5*rho_L*v_L*v_L ;
    E_R = (P_R + gamma_R*p_inf_R)/(gamma_R - 1.0) + 0.5*rho_R*u_R*u_R + 0.5*rho_R*v_R*v_R ;

    S_L = PetscMin((un_R - c_R), (un_L - c_L)) ; 
    S_R = PetscMax((un_L + c_L), (un_R + c_R)) ;
    S_star = ( P_R - P_L + rho_L*un_L*(S_L-un_L) - rho_R*un_R*(S_R - un_R) )/(rho_L*(S_L - un_L) - rho_R*(S_R - un_R)) ;

    // Now compute the left right and starred fluxes for HLLC.
    
    FL[5] = un_L ; FR[5] = un_R ;

    UL_star[0] = rho_L*(S_L - un_L)/(S_L - S_star) ;  
    un = S_star ; ut = ut_L ;
    UL_star[1] = UL_star[0]*(un*nx - ut*ny) ;  
    UL_star[2] = UL_star[0]*(un*ny + ut*nx) ;
    UL_star[3] = UL_star[0]*( (E_L/rho_L) + (S_star - un_L)*(S_star + P_L/(rho_L*(S_L - un_L)) ) )	;  
    UL_star[4] = VL[4]*(S_L - un_L)/(S_L - S_star);	

    UR_star[0] = rho_R*(S_R - un_R)/(S_R - S_star) ;
    un = S_star ; ut = ut_R ;
    UR_star[1] = UR_star[0]*(un*nx - ut*ny) ;
    UR_star[2] = UR_star[0]*(un*ny + ut*nx) ;
    UR_star[3] = UR_star[0]*( (E_R/rho_R) + (S_star - un_R)*(S_star + P_R/(rho_R*(S_R - un_R)) ) ) ;
    UR_star[4] = VR[4]*(S_R - un_R)/(S_R - S_star) ;

    for(i = 0 ; i < nVar ; i++) 
        FL_star[i] = FL[i] + S_L*(UL_star[i] - QL[i]) ; 
    
    FL_star[5] = un_L + S_L*( ((S_L - un_L)/(S_L - S_star)) - 1.0 ) ;
    
    for(i = 0 ; i < nVar ; i++) 
        FR_star[i] = FR[i] + S_R*(UR_star[i] - QR[i]) ; 
    
    FR_star[5] = un_R + S_R*( ((S_R - un_R)/(S_R - S_star)) - 1.0 ) ;
    
    if( S_L > 0.0 ) {
        for(i = 0 ; i < 6; i++) {
            Flux[i] = FL[i] ; 
        }
    } 
    
    else if((S_star >= 0.0) && (S_L < 0.0)) {
        for(i = 0 ; i < 6; i++) { 
            Flux[i] = FL_star[i];
        }
    } 
    
    else if((S_star < 0.0) && (S_R >= 0.0)) {
        
        for(i = 0 ; i < 6; i++) {
            Flux[i] = FR_star[i];
        }
    } 
    
    else if(S_R < 0.0) {
        for(i = 0 ; i < 6; i++) {
            Flux[i] = FR[i]; 
        }
    }
    
    return s_max; 
}

//----------------------------------------------------------------------------
// Rotated HLLC Riemann solver 
//----------------------------------------------------------------------------

PetscReal rotHLLCRiemannSolver(const PetscReal* VL, const PetscReal* VR,
                         const PetscReal nx, const PetscReal ny,
                         const PetscReal x,  const PetscReal y,
                         PetscReal* Flux) {
    
    
    PetscReal flux1[6];
    PetscReal flux2[6];
    
    PetscReal smax1, smax2;

    PetscInt i ;
    PetscReal alpha1, alpha2, n1x, n1y, n2x, n2y, u_L, u_R, v_L, v_R, du, dv, dq ;

    u_L = VL[1]; v_L = VL[2]; 
    u_R = VR[1]; v_R = VR[2]; 
    du = u_R - u_L ; dv = v_R - v_L ; 
    dq = PetscSqrtReal(du*du + dv*dv) ;
    
    if(dq < 1.0E-10) { n1x = nx ; n1y = ny ; }
    else { n1x = du/dq ; n1y = dv/dq ; }

    alpha1 = (n1x*nx + n1y*ny) ;
    if(alpha1 < 0) { n1x = -n1x ; n1y = -n1y ; alpha1 = -alpha1 ; }
    n2x = -n1y ; n2y = n1x ;
    alpha2 = (n2x*nx + n2y*ny) ; 
    if(alpha2 < 0) { n2x = -n2x ; n2y = -n2y ; alpha2 = -alpha2 ; }

    smax1 = HLLCRiemannSolver(VL, VR, n1x, n1y, x, y, flux1);
    smax2 = HLLCRiemannSolver(VL, VR, n2x, n2y, x, y, flux2);

    for(i = 0 ; i < nVar+1; i++) 
        Flux[i] = alpha1*flux1[i] + alpha2*flux2[i];
    
    
    return PetscMax(smax1, smax2);
                                  
}

//----------------------------------------------------------------------------
// Viscous Riemann Solver (Does average of the two fluxes)  
//----------------------------------------------------------------------------

PetscReal ViscLLFRiemannSolverPrim(PetscReal y, const PetscReal* VL, PetscReal grad_VL[nVar][DIM],
                               const PetscReal* VR, PetscReal grad_VR[nVar][DIM], 
                               const PetscReal nx, const PetscReal ny,
                               PetscReal* Flux) {
    
    PetscReal FL[6], FR[6]; 
    PetscInt c; 
    
    PetscReal s_max_l = PDEViscFluxPrim(y, VL, grad_VL, nx, ny, FL);
    PetscReal s_max_r = PDEViscFluxPrim(y, VR, grad_VR, nx, ny, FR);
    
    PetscReal s_max = PetscMax(s_max_l, s_max_r);
    
    for (c = 0; c < 6; ++c) 
        Flux[c] = 0.5*(FR[c] + FL[c]);
    
    return s_max; 
}

