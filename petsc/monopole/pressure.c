/*
 * initial_condition.c
 *      Author: sunder
 */ 
#include "hype.h" 

//----------------------------------------------------------------------------
// Point source as a function of time  
//----------------------------------------------------------------------------

PetscReal Source(PetscReal  t){

    PetscReal   f0 = 100.0;      /*dominant frequency*/
    PetscReal   t0 = 4.0/ f0;    /*source time shift*/
    
    return -2.0*(t - t0)*f0*f0*PetscExpReal( -1.0*f0*f0*(t - t0)*(t - t0));
}


//----------------------------------------------------------------------------------------
// Pressure and its time derivative at point x and time t emitted by a monopole at (0,0,0)
//----------------------------------------------------------------------------------------

Field Pressure(PetscReal x, PetscReal y, PetscReal z, PetscReal t){

    Field Q0;

    PetscReal r = PetscSqrtReal(x*x + y*y + z*z);
    
    if (r != 0.0)
        Q0.comp[0] =  Source(t - r/c0)/(4.0*PETSC_PI*r);
	 else 
         Q0.comp[0] = 0.0;

    return Q0;
}

//----------------------------------------------------------------------------
// Compute exact pressure at every cell at a given time
//----------------------------------------------------------------------------
PetscErrorCode ComputePressureExact(Vec U, DM da, PetscReal t, AppCtx Ctx) {

    PetscErrorCode ierr;

    DM          coordDA;
    Vec         coordinates;
    DMDACoor2d  **coords;
    Field   **u;
    PetscInt    xs, ys, xm, ym, i, j, c;
    PetscReal zc, rc;
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
            
             Q0 = Pressure(rc, 0.0, zc, t);

            for (c = 0; c < nVar; ++c)
                u[j][i].comp[c] = Q0.comp[c];
            
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);

    return ierr; 
}
