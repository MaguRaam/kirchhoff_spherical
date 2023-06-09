/*
 * initialize_solution.c
 *      Author: sunder
 */ 

#include "hype.h" 

//----------------------------------------------------------------------------
// Initialize the solution with the given initial condition  
//----------------------------------------------------------------------------

PetscErrorCode InitializeSolution(Vec U, DM da, AppCtx Ctx) {

    PetscErrorCode ierr;

    DM          coordDA;
    Vec         coordinates;
    DMDACoor2d  **coords;
    Field   **u;
    PetscInt    xs, ys, xm, ym, i, j, c, l , m;
    PetscReal integral[nVar]; 
    PetscReal xc, yc, xGP, yGP;
    PetscReal h = Ctx.h; 
    PetscReal Q0[nVar]; 

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(da, &coordDA);CHKERRQ(ierr);
    ierr = DMGetCoordinates(da, &coordinates);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);CHKERRQ(ierr);

    // Use five point gauss quadrature

    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {
            
            // Get coordinates of center of the cell 
            
            xc = coords[j][i].x; 
            yc = coords[j][i].y;
            
            for (c = 0; c < nVar; ++c)
                integral[c] = 0.0;
            
            for(l = 0; l < N_gp5; ++l) {
                for (m = 0; m < N_gp5; ++m) {

                    xGP = xc + h*x_gp5[l];
                    yGP = yc + h*x_gp5[m];
                    
                    InitialCondition(xGP,yGP,Q0);
                    
                    for (c = 0; c < nVar; ++c) 
                        integral[c] += w_gp5[l]*w_gp5[m]*Q0[c];
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
