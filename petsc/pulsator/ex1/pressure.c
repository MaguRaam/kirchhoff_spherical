/*
 * initial_condition.c
 *      Author: sunder
 */ 
#include "hype.h" 


//----------------------------------------------------------------------------------------
// Pressure emitted by rigid pulsator
//----------------------------------------------------------------------------------------

Field Pressure(PetscReal x, PetscReal y, PetscReal z, PetscReal t){

    Field Q0;

    PetscReal r = PetscSqrtReal(x*x + y*y + z*z);
    
    PetscReal rho0 = 1000.0;             // density of water medium
    PetscReal c = 1437.0;                // speed of sound in water
    PetscReal U0 = 20.0;
    PetscReal lambda = 2.0e-3;           // wavelength
    PetscReal k = (2.0 * M_PI) / lambda; // wavenumber
    PetscReal omega = 4.7837e6;          // frequency of pulsator
    PetscReal R0 = 2.0e-5;               // radius of sphere

    PetscReal cosX0 = (k * R0) / sqrt(1.0 + k * k * R0 * R0);
    PetscReal X0 = acos(cosX0);

    Q0.comp[0] = rho0 * c * U0 * (R0 / r) * cosX0 * exp(-k * (r - R0) + X0) * cos(omega * t);

    return Q0;
}

//----------------------------------------------------------------------------
// Compute cell average of pressure at every cell at a given time
//----------------------------------------------------------------------------
PetscErrorCode ComputePressureAverage(Vec U, DM da, PetscReal t, AppCtx Ctx) {

    PetscErrorCode ierr;

    DM          coordDA;
    Vec         coordinates;
    DMDACoor2d  **coords;
    Field   **u;
    PetscInt    xs, ys, xm, ym, i, j, c, l, m;
    PetscReal zc, rc, zGP, rGP;
    PetscReal h = Ctx.h;
    PetscReal integral[nVar];
    Field     Q0;

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(da, &coordDA);CHKERRQ(ierr);
    ierr = DMGetCoordinates(da, &coordinates);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);CHKERRQ(ierr);

    // Use five point gauss quadrature

    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {
            
            // Get coordinates of center of the cell 
            
            zc = coords[j][i].x; 
            rc = coords[j][i].y;
            
            for (c = 0; c < nVar; ++c)
                integral[c] = 0.0;

             for(l = 0; l < N_gp5; ++l) {
                for (m = 0; m < N_gp5; ++m) {

                    zGP = zc + h*x_gp5[l];
                    rGP = rc + h*x_gp5[m];
                    
                    Q0 = Pressure(0.0, rGP, zGP, t);
                    
                    for (c = 0; c < nVar; ++c) 
                        integral[c] += w_gp5[l]*w_gp5[m]*Q0.comp[c];
                }
            }

            for (c = 0; c < nVar; ++c)
                u[j][i].comp[c] = integral[c];
            
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);

    return ierr; 
}
