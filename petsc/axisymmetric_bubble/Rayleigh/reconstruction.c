/*
 * reconstruction.c
 *      Author: sunder
 */ 
#include "hype.h"

//----------------------------------------------------------------------------
// Value of Nth Order basis functions 
//----------------------------------------------------------------------------

PetscReal basis(PetscReal x, PetscReal y, PetscInt n) {
    
    switch (n) {
        case 0:
            return 1.0;
            break; 
        case 1:
            return x; 
            break;
        case 2:
            return y;
            break;
        case 3:
            return x*x - 1./12.;
            break;
        case 4:
            return y*y - 1./12.;
            break;
        case 5:
            return x*y;
            break; 
        case 6:
            return x*(x*x - 3./20.);
            break;
        case 7:
            return y*(y*y - 3./20.);
            break; 
        case 8:
            return y*(x*x - 1./12.);
            break; 
        case 9:
            return x*(y*y - 1./12.);
            break;
        default:
            return 0.0;
    }
}

//----------------------------------------------------------------------------
// Gradients of Nth Order basis functions 
//----------------------------------------------------------------------------

void basis_grad(PetscReal x, PetscReal y, PetscInt n, PetscReal* grad_x, PetscReal* grad_y) {
    
    switch (n) {
        case 0:
            *grad_x = 0.0;
            *grad_y = 0.0; 
            break; 
        case 1:
            *grad_x = 1.0;
            *grad_y = 0.0; 
            break;
        case 2:
            *grad_x = 0.0;
            *grad_y = 1.0; 
            break;
        case 3:
            *grad_x = 2.0*x;
            *grad_y = 0.0; 
            break;
        case 4:
            *grad_x = 0.0;
            *grad_y = 2.0*y; 
            break;
        case 5:
            *grad_x = y;
            *grad_y = x; 
            break; 
        case 6:
            *grad_x = 3.0*x*x - 3./20.;
            *grad_y = 0.0; ;
            break;
        case 7:
            *grad_x = 0.0;
            *grad_y = 3.0*y*y - 3./20.;
            break; 
        case 8:
            *grad_x = 2.0*x*y;
            *grad_y = (x*x - 1./12.); 
            break; 
        case 9:
            *grad_x = (y*y - 1./12.);
            *grad_y = 2.0*x*y;
            break;
        default:
            *grad_x = 0.0;
            *grad_y = 0.0; 
    }
}


PetscReal evaluate_polynomial(const PetscReal x, const PetscReal y, const PetscReal coeffs[]) {
  
  return coeffs[0] + coeffs[1]*x + coeffs[2]*y + coeffs[3]*(x*x - 1./12.) + coeffs[4]*(y*y - 1./12.) +
         coeffs[5]*(x*y) + coeffs[6]*(x*(x*x - 3./20.)) + coeffs[7]*(y*(y*y - 3./20.)) + coeffs[8]*(y*(x*x - 1./12.)) + coeffs[9]*(x*(y*y - 1./12.));
}   


void evaluate_grad(const PetscReal coeffs[], PetscReal x, PetscReal y, const PetscReal h, 
                                PetscReal* grad_x, PetscReal* grad_y) {

  *grad_x = coeffs[1] + coeffs[3]*(2.0*x) + coeffs[5]*y + coeffs[6]*(3.0*x*x - 3./20.) + coeffs[8]*(2.0*x*y) + coeffs[9]*(y*y - 1./12.);
  *grad_y = coeffs[2] + coeffs[4]*(2.0*y) + coeffs[5]*x + coeffs[7]*(3.0*y*y - 3./20.) + coeffs[8]*((x*x - 1./12.)) + coeffs[9]*(2.0*x*y);
}






//----------------------------------------------------------------------------
// Various slope limiters 
//----------------------------------------------------------------------------

PetscReal minmod(PetscReal a, PetscReal b) {
    if (a*b < 0.0) {
        return 0.0;
    }
    
    else {
        if (PetscAbsReal(a) < PetscAbsReal(b) )
            return a;
        else 
            return b; 
    }
}

//----------------------------------------------------------------------------
// 2D 4th order WENO reconstruction 
//----------------------------------------------------------------------------

PetscReal pow4(PetscReal a) {
    PetscReal a2 = a*a; 
    return a2*a2; 
}

void weno(const PetscReal U_x[], const PetscReal U_y[], const PetscReal U_xy[], PetscReal u_coeffs[], PetscReal coeffs[]) {
    
    PetscReal gammaHi = 0.85; 
    PetscReal gammaLo = 0.85; 
    PetscReal u_0 = U_xy[0]; 
    PetscReal u_ip1 = U_x[3]; PetscReal u_jp1 = U_y[3]; PetscReal u_ip1jp1 = U_xy[1];
    PetscReal u_im1 = U_x[1]; PetscReal u_jm1 = U_y[1]; PetscReal u_ip1jm1 = U_xy[2];
    PetscReal u_ip2 = U_x[4]; PetscReal u_jp2 = U_y[4]; PetscReal u_im1jp1 = U_xy[3];
    PetscReal u_im2 = U_x[0]; PetscReal u_jm2 = U_y[0]; PetscReal u_im1jm1 = U_xy[4];
    PetscReal u_xR4, u_yR4, u_xxR4, u_yyR4, u_xyR4, u_xxxR4, u_yyyR4, u_xxyR4, u_xyyR4;
    PetscReal u_xR3[5]; PetscReal u_yR3[5]; PetscReal u_xxR3[5]; PetscReal u_yyR3[5]; PetscReal u_xyR3[5];
    PetscReal IS_R4; PetscReal gamma_R4; PetscReal w_R4; 
    PetscReal IS_R3[5]; PetscReal gamma_R3[5];  PetscReal w_R3[5];  PetscReal sum = 0.0; 
    PetscReal wt_ratio; 
    PetscInt i; 

    gamma_R4 = gammaHi;

    gamma_R3[0] = (1.0 - gammaHi)*gammaLo; 

    for (i = 1 ; i < 5; ++i)
        gamma_R3[i] = 0.25*(1-gammaHi)*(1-gammaLo);

    // Fourth order stencil 

    u_xR4   = r41_60*(-u_im1 + u_ip1) + r11_120*(u_im2 - u_ip2);
    u_yR4   = r41_60*(-u_jm1 + u_jp1) + r11_120*(u_jm2 - u_jp2);
    u_xxR4  = -u_0 + 0.5*(u_im1 + u_ip1);
    u_yyR4  = -u_0 + 0.5*(u_jm1 + u_jp1);
    u_xyR4  = 0.25*(u_im1jm1 - u_im1jp1 - u_ip1jm1 + u_ip1jp1); 
    u_xxxR4 = r1_6*(u_im1 - u_ip1) + r1_12*(u_ip2 - u_im2); 
    u_yyyR4 = r1_6*(u_jm1 - u_jp1) + r1_12*(u_jp2 - u_jm2);
    u_xxyR4 = 0.25*(-u_im1jm1 + u_im1jp1 - u_ip1jm1 + u_ip1jp1) + 0.5*(u_jm1 - u_jp1);
    u_xyyR4 = 0.5*(u_im1 - u_ip1) + 0.25*(u_ip1jp1 - u_im1jm1 - u_im1jp1  + u_ip1jm1); 

    u_coeffs[0] = u_0;
    u_coeffs[1] = u_xR4;
    u_coeffs[2] = u_yR4;
    u_coeffs[3] = u_xxR4;
    u_coeffs[4] = u_yyR4;
    u_coeffs[5] = u_xyR4;
    u_coeffs[6] = u_xxxR4;
    u_coeffs[7] = u_yyyR4;
    u_coeffs[8] = u_xxyR4;
    u_coeffs[9] = u_xyyR4;

    IS_R4 = (u_xR4 + 0.1*u_xxxR4)*(u_xR4 + 0.1*u_xxxR4) + r13_3*(u_xxR4*u_xxR4 + u_yyR4*u_yyR4) + 
            (u_yR4 + 0.1*u_yyyR4)*(u_yR4 + 0.1*u_yyyR4) + 39.05*(u_xxxR4*u_xxxR4 + u_yyyR4*u_yyyR4) + 
            r7_6*u_xyR4*u_xyR4 + 4.7*(u_xxyR4*u_xxyR4 + u_xyyR4*u_xyyR4);
            
    w_R4 = gamma_R4/(pow4(IS_R4 + small_num)); 
    sum = w_R4; 

    // Stencil 1 (centered stencil)

    u_xR3[0]  = 0.5*(-u_im1 + u_ip1);
    u_yR3[0]  = 0.5*(-u_jm1 + u_jp1);
    u_xxR3[0] = -u_0 + 0.5*(u_im1 + u_ip1);
    u_yyR3[0] = -u_0 + 0.5*(u_jm1 + u_jp1);
    u_xyR3[0] = 0.25*(u_im1jm1 - u_im1jp1 - u_ip1jm1 + u_ip1jp1); 

    // Stencil 2 (left biased stencil)

    u_xR3[1]  = 1.5*u_0 - 2.0*u_im1 + 0.5*u_im2;
    u_yR3[1]  = 0.5*(-u_jm1 + u_jp1); 
    u_xxR3[1] = 0.5*u_0 - u_im1 + 0.5*u_im2; 
    u_yyR3[1] = -u_0 + 0.5*u_jm1 + 0.5*u_jp1;
    u_xyR3[1] = 0.5*(u_im1jm1 - u_im1jp1 - u_jm1 + u_jp1);

    // Stencil 3 (right biased stencil)

    u_xR3[2]  = -1.5*u_0 + 2.0*u_ip1 - 0.5*u_ip2;
    u_yR3[2]  = 0.5*(-u_jm1 + u_jp1);
    u_xxR3[2] = 0.5*u_0 - u_ip1 + 0.5*u_ip2; 
    u_yyR3[2] = -u_0 + 0.5*u_jm1 + 0.5*u_jp1;
    u_xyR3[2] = 0.5*(-u_ip1jm1 + u_ip1jp1 + u_jm1 - u_jp1);

    // Stencil 4 (bottom biased stencil)

    u_xR3[3]  = 0.5*(-u_im1 + u_ip1);
    u_yR3[3]  = 1.5*u_0 - 2.0*u_jm1 + 0.5*u_jm2; 
    u_xxR3[3] = -u_0 + 0.5*u_im1 + 0.5*u_ip1;
    u_yyR3[3] = 0.5*u_0 - u_jm1 + 0.5*u_jm2;
    u_xyR3[3] = 0.5*(-u_im1 + u_im1jm1 + u_ip1 - u_ip1jm1); 

    // Stencil 5 (top biased stencil)

    u_xR3[4]  = 0.5*(-u_im1 + u_ip1);
    u_yR3[4]  = -1.5*u_0 + 2.0*u_jp1 - 0.5*u_jp2;
    u_xxR3[4] = -u_0 + 0.5*u_im1 + 0.5*u_ip1;
    u_yyR3[4] = 0.5*u_0 - u_jp1 + 0.5*u_jp2; 
    u_xyR3[4] = 0.5*(u_im1 - u_im1jp1 - u_ip1 + u_ip1jp1);

    // Find the smoothness indicators 

    for (i = 0; i < 5; ++i) {
        IS_R3[i] = u_xR3[i]*u_xR3[i] + u_yR3[i]*u_yR3[i] + r13_3*(u_xxR3[i]*u_xxR3[i] + u_yyR3[i]*u_yyR3[i]) + r7_6*u_xyR3[i]*u_xyR3[i];
        w_R3[i] = gamma_R3[i]/(pow4(IS_R3[i] + small_num)); 
        sum += w_R3[i]; 
    }

    // Normalize the weights 

    w_R4 = w_R4/sum; 

    for (i = 0; i < 5; ++i)
        w_R3[i] = w_R3[i]/sum; 

    wt_ratio = w_R4/gamma_R4;

    coeffs[0] = u_0; 

    coeffs[1] = wt_ratio*(u_xR4 - gamma_R3[0]*u_xR3[0] - gamma_R3[1]*u_xR3[1] - gamma_R3[2]*u_xR3[2] - gamma_R3[3]*u_xR3[3] - gamma_R3[4]*u_xR3[4])
                                + w_R3[0]*u_xR3[0] + w_R3[1]*u_xR3[1] + w_R3[2]*u_xR3[2] + w_R3[3]*u_xR3[3] + w_R3[4]*u_xR3[4]; 

    coeffs[2] = wt_ratio*(u_yR4 - gamma_R3[0]*u_yR3[0] - gamma_R3[1]*u_yR3[1] - gamma_R3[2]*u_yR3[2] - gamma_R3[3]*u_yR3[3] - gamma_R3[4]*u_yR3[4])
                                + w_R3[0]*u_yR3[0] + w_R3[1]*u_yR3[1] + w_R3[2]*u_yR3[2] + w_R3[3]*u_yR3[3] + w_R3[4]*u_yR3[4];
                                
    coeffs[3] = wt_ratio*(u_xxR4 - gamma_R3[0]*u_xxR3[0] - gamma_R3[1]*u_xxR3[1] - gamma_R3[2]*u_xxR3[2] - gamma_R3[3]*u_xxR3[3] - gamma_R3[4]*u_xxR3[4])
                                + w_R3[0]*u_xxR3[0] + w_R3[1]*u_xxR3[1] + w_R3[2]*u_xxR3[2] + w_R3[3]*u_xxR3[3] + w_R3[4]*u_xxR3[4];
                                
    coeffs[4] = wt_ratio*(u_yyR4 - gamma_R3[0]*u_yyR3[0] - gamma_R3[1]*u_yyR3[1] - gamma_R3[2]*u_yyR3[2] - gamma_R3[3]*u_yyR3[3] - gamma_R3[4]*u_yyR3[4])
                                + w_R3[0]*u_yyR3[0] + w_R3[1]*u_yyR3[1] + w_R3[2]*u_yyR3[2] + w_R3[3]*u_yyR3[3] + w_R3[4]*u_yyR3[4];
                                
    coeffs[5] = wt_ratio*(u_xyR4 - gamma_R3[0]*u_xyR3[0] - gamma_R3[1]*u_xyR3[1] - gamma_R3[2]*u_xyR3[2] - gamma_R3[3]*u_xyR3[3] - gamma_R3[4]*u_xyR3[4])
                                + w_R3[0]*u_xyR3[0] + w_R3[1]*u_xyR3[1] + w_R3[2]*u_xyR3[2] + w_R3[3]*u_xyR3[3] + w_R3[4]*u_xyR3[4];
                                
    coeffs[6] = wt_ratio*u_xxxR4;

    coeffs[7] = wt_ratio*u_yyyR4;

    coeffs[8] = wt_ratio*u_xxyR4;

    coeffs[9] = wt_ratio*u_xyyR4;
}

