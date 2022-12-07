/*
 * monitor_function.c
 *      Author: sunder
 */ 

#include "hype.h"  

//----------------------------------------------------------------------------
// Write pressure data in vts format
//----------------------------------------------------------------------------
PetscErrorCode WriteVtk(DM da, Vec U, PetscInt step)
{
    PetscErrorCode ierr; 

    char filename[30];
    PetscViewer viewer;

    sprintf(filename, "plot/sol-%08d.vts", step); // 8 is the padding level, increase it for longer simulations
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
    ierr = DMView(da, viewer);
    VecView(U, viewer);

    ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
    
    return ierr; 
}