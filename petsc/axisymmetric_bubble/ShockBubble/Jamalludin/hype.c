/*
 * hype.c
 *      Author: sunder
 */ 

static char help[] = "4th Order 2D code for solving Multiphase Euler equations using PETSc.\n\n";

#include "hype.h" 

//----------------------------------------------------------------------------
// Main function of the code 
//----------------------------------------------------------------------------

int main(int argc,char **argv) {
    
    // --------------------------------------------
    // Initialize MPI 
    //---------------------------------------------

    PetscErrorCode ierr;                    /* For catching PETSc errors */ 
    PetscLogDouble start_time, end_time;    /* For logging the time values */

    ierr = PetscInitialize(&argc, &argv, (char*)0, help);CHKERRQ(ierr);
    ierr =  PetscTime(&start_time);CHKERRQ(ierr); 

    // --------------------------------------------
    // Set important user defined parameters  
    //---------------------------------------------

    AppCtx Ctx; 

    Ctx.x_min           = -6.6e-4;
    Ctx.x_max           =  3.4e-4;
    Ctx.y_min           =  0.0;
    Ctx.y_max           =  5.e-4;
    Ctx.N_x             =  392;
    Ctx.N_y             =  196;
    Ctx.CFL             =  0.9;
    Ctx.InitialStep     =  0;
    Ctx.InitialTime     =  0.0;
    Ctx.FinalTime       =  1.e-6;
    Ctx.WriteInterval   =  10;      
    Ctx.RestartInterval =  500;
    Ctx.left_boundary   =  transmissive;
    Ctx.right_boundary  =  transmissive;
    Ctx.bottom_boundary =  transmissive;
    Ctx.top_boundary    =  transmissive;
    Ctx.Restart         =  PETSC_FALSE;
    Ctx.h               =  (Ctx.x_max - Ctx.x_min)/(PetscReal)(Ctx.N_x);

    //set kirchhoff parameters:

    Ctx.KirchhoffWriteInterval = 10;
    Ctx.surface.Rk      =  5.0*R0;
    Ctx.surface.Nk      =  90;
    
    //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // No need to change anything beyond this point 
    //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // --------------------------------------------
    // Data members  
    //---------------------------------------------

    Vec U;                           // Solution Vector (Conserved variables)
    Vec RHS;                         // RHS vector to update the solution 
    DM da;                           // Grid object 
    PetscInt time_steps;             // No. of time steps 
    TS ts;                           // Time stepping object 
    PetscMPIInt MyPID;               // Rank of the current processor 
    PetscMPIInt numProcs;            // Size of the communicator

    // --------------------------------------------
    // Obtain the rank of the process and size of 
    // the communicator 
    //---------------------------------------------

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcs);CHKERRQ(ierr); 
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Code running with %d processes\n", numProcs);CHKERRQ(ierr); 

    // --------------------------------------------
    // Initialize the grid and set field names
    //---------------------------------------------

    DMBoundaryType x_boundary;
    DMBoundaryType y_boundary;

    if (Ctx.left_boundary == periodic || Ctx.right_boundary == periodic)
        x_boundary = DM_BOUNDARY_PERIODIC;
    else
        x_boundary = DM_BOUNDARY_GHOSTED; 

    if (Ctx.bottom_boundary == periodic || Ctx.top_boundary == periodic)
        y_boundary = DM_BOUNDARY_PERIODIC;
    else
        y_boundary = DM_BOUNDARY_GHOSTED; 

    ierr = DMDACreate2d(PETSC_COMM_WORLD, // Global communicator      
                        x_boundary,       // Boundary conditions in x-direction 
                        y_boundary,       // Boundary conditions in y-direction
                        DMDA_STENCIL_BOX, // Stencil type (other is star type)
                        Ctx.N_x,          // No. of cells in x-direction 
                        Ctx.N_y,          // No. of cells in y-direction
                        PETSC_DECIDE,     // Domain decomposition in x-direction 
                        PETSC_DECIDE,     // Domain decomposition in y-direction
                        nVar,             // No. of dofs per cell 
                        3,                // Width of the stencil
                        NULL,
                        NULL,
                        &da);CHKERRQ(ierr); // da object 

    ierr = DMSetUp(da);CHKERRQ(ierr);

    // Now create various global vectors 

    ierr = DMCreateGlobalVector(da, &U);CHKERRQ(ierr);
    ierr = VecDuplicate(U,&Ctx.W);CHKERRQ(ierr);
    ierr = VecDuplicate(U,&RHS);CHKERRQ(ierr);

    // Set coordinates of cell centers 

    ierr = DMDASetUniformCoordinates(da,
                                        Ctx.x_min + 0.5*Ctx.h, Ctx.x_max - 0.5*Ctx.h,
                                        Ctx.y_min + 0.5*Ctx.h, Ctx.y_max - 0.5*Ctx.h,
                                        0.0,0.0);CHKERRQ(ierr);
   
   
    // Kirchhoff:
    ierr = ComputeNofQuadPoints(da, &Ctx); 
    ierr = AllocateKirchhoffData(da, &Ctx);

    // Set names of the fields
                                        
    ierr = PetscObjectSetName((PetscObject)U,"cons");CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)Ctx.W,"sol");CHKERRQ(ierr);

    ierr = DMDASetFieldName(da,0,"density");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,1,"z-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,2,"r-velocity");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,3,"pressure");CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,4,"alpha-1");CHKERRQ(ierr);
    
    // --------------------------------------------
    // Allocate memory for boundary values and 
    // upwind fluxes
    //---------------------------------------------

    PetscInt xs,ys,xm,ym,q,k;
    PetscReal value, grad_x, grad_y;

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    Ctx.u_bnd        = allocate5d(ym+2, xm+2, nVar, 4, N_gp2);  // 4->number of faces in a cell
    Ctx.u_bnd_grad   = allocate6d(ym+2, xm+2, nVar, 4, N_gp2, DIM); // 4->number of faces, 2->number of quadraure points
    Ctx.F            = allocate3d(ym, xm+1, nVar+1);
    Ctx.G            = allocate3d(ym+1, xm, nVar+1);
    Ctx.phiFace      = allocate3d(4, N_gp2, nDOF);     // 4->number of faces in a cell
    Ctx.phiNode      = allocate2d(N_node, nDOF);
    Ctx.phiVol       = allocate2d(N_gp2d, nDOF);
    Ctx.gradphiVol_x = allocate2d(N_gp2d, nDOF);
    Ctx.gradphiVol_y = allocate2d(N_gp2d, nDOF);
    Ctx.gradphiFace_x = allocate3d(4, N_gp2, nDOF);
    Ctx.gradphiFace_y = allocate3d(4, N_gp2, nDOF);
    
    // Find the value of basis functions on the face quadrature points  
    
    for (q = 0; q < N_gp2; ++q) {
        for (k = 0; k < nDOF; ++k) {
            
            // Left face 
            value = basis(-0.5, x_gp2[q], k);
            basis_grad(-0.5, x_gp2[q], k, &grad_x, &grad_y);
            set_element_3d(Ctx.phiFace, 0, q, k, value);
            set_element_3d(Ctx.gradphiFace_x, 0, q, k, grad_x);
            set_element_3d(Ctx.gradphiFace_y, 0, q, k, grad_y);
            
            // Right face 
            value = basis( 0.5, x_gp2[q], k);
            basis_grad( 0.5, x_gp2[q], k, &grad_x, &grad_y);
            set_element_3d(Ctx.phiFace, 1, q, k, value);
            set_element_3d(Ctx.gradphiFace_x, 1, q, k, grad_x);
            set_element_3d(Ctx.gradphiFace_y, 1, q, k, grad_y);
            
            // Bottom face 
            value = basis(x_gp2[q], -0.5, k);
            basis_grad(x_gp2[q], -0.5, k, &grad_x, &grad_y);
            set_element_3d(Ctx.phiFace, 2, q, k, value);
            set_element_3d(Ctx.gradphiFace_x, 2, q, k, grad_x);
            set_element_3d(Ctx.gradphiFace_y, 2, q, k, grad_y);
            
            // Top face 
            value = basis(x_gp2[q], 0.5, k);
            basis_grad(x_gp2[q], 0.5, k, &grad_x, &grad_y);
            set_element_3d(Ctx.phiFace, 3, q, k, value);
            set_element_3d(Ctx.gradphiFace_x, 3, q, k, grad_x);
            set_element_3d(Ctx.gradphiFace_y, 3, q, k, grad_y);

            
        }
    }
    
    // Find the value of basis functions on interior nodes 
    
    for (q = 0; q < N_node; ++q) {
        for (k = 0; k < nDOF; ++k) {
            set_element_2d(Ctx.phiNode, q, k, basis(x_node[q], y_node[q], k));
        }
    }

    // Find the value of basis functions and gradients on interior quadrature points

    for (q = 0; q < N_gp2d; ++q) {
        for (k = 0; k < nDOF; ++k) {
            set_element_2d(Ctx.phiVol, q, k, basis(x_gp2d[q], y_gp2d[q], k));
            basis_grad(x_gp2d[q], y_gp2d[q], k, &grad_x, &grad_y);
            set_element_2d(Ctx.gradphiVol_x, q, k, grad_x);
            set_element_2d(Ctx.gradphiVol_y, q, k, grad_y);
        }
    }

    ierr = DMCreateLocalVector(da,&Ctx.localU);CHKERRQ(ierr);

    // --------------------------------------------
    // Initialize the solution (either with initial
    // condition or restart file)
    //---------------------------------------------

    if (Ctx.Restart) {
        
        // Initialize by reading the restart file 
        
        PetscViewer    viewer_binary;
        ierr = PetscPrintf(PETSC_COMM_WORLD,"reading vector in binary from restart1.bin ...\n");CHKERRQ(ierr);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"restart1.bin",FILE_MODE_READ,&viewer_binary);CHKERRQ(ierr);
        ierr = VecLoad(U,viewer_binary);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer_binary);CHKERRQ(ierr);
    }

    else {
        
        // Initialize by initial condition 
        ierr = InitializeSolution(U, da, Ctx);CHKERRQ(ierr);
    }

    // --------------------------------------------
    // Advance solution in time   
    //---------------------------------------------
    
    ierr = TSCreate(PETSC_COMM_SELF, &ts);CHKERRQ(ierr);              
    ierr = TSSetProblemType(ts,TS_NONLINEAR);CHKERRQ(ierr);           
    ierr = TSSetDM(ts,da);CHKERRQ(ierr);                              
    
    ierr = RHSFunction(ts, Ctx.InitialTime, U, RHS, &Ctx);CHKERRQ(ierr);
    ierr = TSSetRHSFunction(ts,NULL,RHSFunction, &Ctx);CHKERRQ(ierr);
    
    ierr = TSSetStepNumber(ts,Ctx.InitialStep);
    ierr = TSSetTime(ts, Ctx.InitialTime);
    ierr = TSSetTimeStep(ts, Ctx.dt);CHKERRQ(ierr);                     
    ierr = TSMonitorSet(ts,MonitorFunction,&Ctx,NULL);CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts, Ctx.FinalTime);CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP);CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts);CHKERRQ(ierr); 
    ierr = TSSetType(ts, TSSSP);CHKERRQ(ierr); 
    //ierr = TSSSPSetType(ts, TSSSPRK104);CHKERRQ(ierr);
    ierr = TSSSPSetType(ts, TSSSPRKS3);CHKERRQ(ierr); 
    ierr = TSSSPSetNumStages(ts,4);CHKERRQ(ierr);; 
    ierr = TSSolve(ts, U);CHKERRQ(ierr);
    ierr = TSGetStepNumber(ts,&time_steps);CHKERRQ(ierr); 
    
    // --------------------------------------------
    // Output solution in vtk format   
    //--------------------------------------------
    
    ierr = ComputePrimitiveVariables(U, Ctx.W, da);CHKERRQ(ierr);
    char filename[20]; 
    sprintf(filename, "sol-%08d.vts", time_steps);
    PetscViewer viewer;  
    ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
    ierr = DMView(da, viewer);CHKERRQ(ierr);
    ierr = VecView(Ctx.W, viewer);CHKERRQ(ierr);
    
    // --------------------------------------------
    // Get the norms of errors (only for periodic
    // test cases)
    //---------------------------------------------

    PetscReal nrm_2, nrm_inf;
    ierr = ErrorNorms(U, da, Ctx, &nrm_2, &nrm_inf);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Norm2 = %.7e, NormMax = %.7e\n", nrm_2, nrm_inf);CHKERRQ(ierr);

    // --------------------------------------------
    // Print the time taken for simulation       
    //---------------------------------------------

    ierr =  PetscTime(&end_time);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time taken =  %g\n",(double)(end_time - start_time));CHKERRQ(ierr);

    // --------------------------------------------
    // Free all the memory, finalize MPI and exit   
    //---------------------------------------------

    ierr = VecDestroy(&U);CHKERRQ(ierr);
    ierr = VecDestroy(&RHS);CHKERRQ(ierr);
    ierr = DMDestroy(&da);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    ierr = TSDestroy(&ts);CHKERRQ(ierr);
    ierr = VecDestroy(&Ctx.localU);CHKERRQ(ierr);
    ierr = VecDestroy(&Ctx.W);CHKERRQ(ierr);
    
    free5d(Ctx.u_bnd);
    free6d(Ctx.u_bnd_grad);
    free3d(Ctx.F);
    free3d(Ctx.G); 
    free3d(Ctx.phiFace);
    free2d(Ctx.phiNode);
    free2d(Ctx.phiVol);
    free2d(Ctx.gradphiVol_x);
    free2d(Ctx.gradphiVol_y);
    free3d(Ctx.gradphiFace_x);
    free3d(Ctx.gradphiFace_y);

    //Kirhhoff:
    ierr = FreeKirchhoffData(&Ctx);

    ierr = PetscFinalize();CHKERRQ(ierr);

    return ierr; 
}
