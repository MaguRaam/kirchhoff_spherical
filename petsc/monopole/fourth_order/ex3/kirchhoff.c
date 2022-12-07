/*
 * kirchhoff.c
 *      Author: Magu
 */ 
#include "hype.h"

//given an aribirtary point x find cell index i:
PetscInt GetCellIndex(PetscReal x, PetscReal xmin, PetscReal h){
   return floor(fabs((x - xmin)/h)); 
}

//----------------------------------------------------------------------------
// Compute no of quadrature points in the current process: 
//----------------------------------------------------------------------------
PetscErrorCode ComputeNofQuadPoints(DM da, AppCtx* Ctx){

    PetscErrorCode  ierr;
    PetscInt        xs, ys, xm, ym, i, j, k;
    PetscReal       thetac, theta1, theta2, xq1, yq1, xq2, yq2;
    PetscReal       Rk = Ctx->surface.Rk;
    PetscInt        Nk = Ctx->surface.Nk;

    // compute dtheta:
    Ctx->surface.dtheta = PETSC_PI / (PetscReal)Nk;
    PetscReal dtheta = Ctx->surface.dtheta;

    // initialize no of quadpoints:
    Ctx->surface.Nk_local = 0;

    // local domain boundaries:
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    // loop over cells on Kirchhoff arc:
    for (k = 0; k < Nk; ++k)
    {
        //compute theta at cell center:
        thetac = (k + 0.5) * dtheta;

        //get theta for quadarture points q1 and q2:
        theta1 = thetac + dtheta*x_gp2[0];
        theta2 = thetac + dtheta*x_gp2[1];

        // coordinates of quadrature point q1:
        xq1 = Rk * PetscCosReal(theta1); yq1 = Rk * PetscSinReal(theta1);

        // coordinates of quadrature point q2:
        xq2 = Rk * PetscCosReal(theta2); yq2 = Rk * PetscSinReal(theta2); 

        // get cell index of quadpoint q1: 
        i = GetCellIndex(xq1, Ctx->x_min, Ctx->h);
        j = GetCellIndex(yq1, Ctx->y_min, Ctx->h);

        //check if the point is in the current process:
        if (i >= xs &&  i < xs+xm && j >= ys &&  j < ys+ym)
            Ctx->surface.Nk_local++;
        
        // get cell index of quadpoint q2: 
        i = GetCellIndex(xq2, Ctx->x_min, Ctx->h);
        j = GetCellIndex(yq2, Ctx->y_min, Ctx->h);

        //check if the point is in the current process:
        if (i >= xs &&  i < xs+xm && j >= ys &&  j < ys+ym)
            Ctx->surface.Nk_local++; 
    }

    return ierr;
}

//----------------------------------------------------------------------------
// Allocate memory for Kirchhoff data 
//----------------------------------------------------------------------------
PetscErrorCode AllocateKirchhoffData(DM da, AppCtx* Ctx){

    PetscErrorCode ierr;
    PetscInt xs, ys, xm, ym, i, j, k;
    PetscInt Nk = Ctx->surface.Nk;
    PetscInt Nk_local = Ctx->surface.Nk_local;
    PetscReal Rk = Ctx->surface.Rk;
    PetscReal dtheta = Ctx->surface.dtheta;
    PetscReal thetac, theta1, theta2, xq1, yq1, xq2, yq2;
    PetscReal distanceq1, distanceq2; 
    PetscInt rank, q;


    //allocate memory for local arary: where size = no of quadpoints in the current process
    Ctx->surface.q_local = (KirchhoffData*)malloc(Nk_local*sizeof(KirchhoffData));
    Ctx->surface.q_point = (Point*)malloc(Nk_local*sizeof(Point));
    Ctx->surface.q_normal = (Point*)malloc(Nk_local*sizeof(Point));
    Ctx->surface.q_cell  = (CellIndex*)malloc(Nk_local*sizeof(CellIndex));     

    //allocate memory in rank 0 for global data:
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr); 
    
    if (rank == 0)
        Ctx->surface.q_global = (KirchhoffData*)malloc(2*Nk*sizeof(KirchhoffData));

    // local domain boundaries:
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    // loop over cells on Kirchhoff arc:
    q = 0;

    for (k = 0; k < Nk; ++k)
    {
        //compute theta at cell center:
        thetac = (k + 0.5) * dtheta;

        //get theta for quadarture points q1 and q2:
        theta1 = thetac + dtheta*x_gp2[0];
        theta2 = thetac + dtheta*x_gp2[1];

        // coordinates of quadrature point q1:
        xq1 = Rk * PetscCosReal(theta1); yq1 = Rk * PetscSinReal(theta1);

        // coordinates of quadrature point q2:
        xq2 = Rk * PetscCosReal(theta2); yq2 = Rk * PetscSinReal(theta2); 
        
        // get cell index of quadpoint q1: 
        i = GetCellIndex(xq1, Ctx->x_min, Ctx->h);
        j = GetCellIndex(yq1, Ctx->y_min, Ctx->h);

        //check if the quad point q1 is in the current process:
        if (i >= xs &&  i < xs+xm && j >= ys &&  j < ys+ym)
        {
            //store quadrature point q1:
            Ctx->surface.q_point[q].x = xq1;
            Ctx->surface.q_point[q].y = yq1;

            //compute distance from center for quadrature point:
            distanceq1 = PetscSqrtReal(xq1*xq1 + yq1*yq1);
            
            //store normal vector for quadrature point q1:
            Ctx->surface.q_normal[q].x = xq1/distanceq1;
            Ctx->surface.q_normal[q].y = yq1/distanceq1;

            //store cell index of quadrature point:
            Ctx->surface.q_cell[q].i = i;
            Ctx->surface.q_cell[q].j = j;

            //store theta:
            Ctx->surface.q_local[q].theta = theta1;

            q++;
        }

        // get cell index of quadpoint q2: 
        i = GetCellIndex(xq2, Ctx->x_min, Ctx->h);
        j = GetCellIndex(yq2, Ctx->y_min, Ctx->h);

        //check if the quad point q2 is in the current process:
        if (i >= xs &&  i < xs+xm && j >= ys &&  j < ys+ym)
        {
            //store quadrature point q2:
            Ctx->surface.q_point[q].x = xq2;
            Ctx->surface.q_point[q].y = yq2;

            //compute distance from center for quadrature point:
            distanceq2 = PetscSqrtReal(xq2*xq2 + yq2*yq2);
            
            //store normal vector for quadrature point q2:
            Ctx->surface.q_normal[q].x = xq2/distanceq2;
            Ctx->surface.q_normal[q].y = yq2/distanceq2;

            //store cell index of quadrature point:
            Ctx->surface.q_cell[q].i = i;
            Ctx->surface.q_cell[q].j = j;

            //store theta:
            Ctx->surface.q_local[q].theta = theta2;

            q++;
        }

    }
    
    return ierr;
}

//----------------------------------------------------------------------------------------------------
// Reconstruct polynomial on a given cell
//----------------------------------------------------------------------------------------------------
void ReconstructOnGivenCell(Field** u, PetscInt i, PetscInt j, PetscReal coeffs[], PetscReal pcoeffs[]){

    PetscReal p_x_loc[s_width], p_y_loc[s_width], p_xy_loc[s_width];

    //get perturbed pressure values from the stencil:
    p_x_loc[0] = u[j][i-2].comp[0]; 
    p_x_loc[1] = u[j][i-1].comp[0]; 
    p_x_loc[2] = u[j][i].comp[0]; 
    p_x_loc[3] = u[j][i+1].comp[0]; 
    p_x_loc[4] = u[j][i+2].comp[0];

    p_y_loc[0] = u[j-2][i].comp[0]; 
    p_y_loc[1] = u[j-1][i].comp[0]; 
    p_y_loc[2] = u[j][i].comp[0]; 
    p_y_loc[3] = u[j+1][i].comp[0]; 
    p_y_loc[4] = u[j+2][i].comp[0];

    p_xy_loc[0] = u[j][i].comp[0];
    p_xy_loc[1] = u[j+1][i+1].comp[0];
    p_xy_loc[2] = u[j-1][i+1].comp[0];
    p_xy_loc[3] = u[j+1][i-1].comp[0];
    p_xy_loc[4] = u[j-1][i-1].comp[0];


    //reconstruct polynomial (get polynomial coefficients)
    weno(p_x_loc, p_y_loc, p_xy_loc, pcoeffs, coeffs);

}

//----------------------------------------------------------------------------
// Write Kirchhoff data 
//----------------------------------------------------------------------------
PetscErrorCode WriteKirchhoffData(Vec U, DM da, AppCtx* Ctx, PetscReal time){

    PetscErrorCode  ierr;
    PetscInt q, i, j, k;
    int nproc, rank;
    PetscReal coeffs[nDOF], pcoeffs[nDOF];
    PetscReal psi, eta, px, py, nx, ny, h = Ctx->h;
    PetscReal xc, yc;
    Field **u;

    // get rank and no of process:
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&nproc);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr); 

    // Scatter global->local to have access to the required ghost values 
    ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, Ctx->localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES,   Ctx->localU);CHKERRQ(ierr);

    // Read the local solution to the array u  
    ierr = DMDAVecGetArrayRead(da, Ctx->localU, &u); CHKERRQ(ierr); 

    // loop over quadarature points in the current process:
    for (q = 0; q < Ctx->surface.Nk_local; ++q){

        // compute cell index containing quad point:
        i = Ctx->surface.q_cell[q].i;
        j = Ctx->surface.q_cell[q].j;

        // reconstruct polynomial on a cell containing quadrature point:
        ReconstructOnGivenCell(u, i, j, coeffs, pcoeffs);

        //compute cell center
        xc = Ctx->x_min + (i + 0.5)*h;
        yc = Ctx->y_min + (j + 0.5)*h;

        //compute reference coordinate of quadrature point:
        psi = (Ctx->surface.q_point[q].x - xc)/h;
        eta = (Ctx->surface.q_point[q].y - yc)/h;

        //compute pressure at quadrature point:
        Ctx->surface.q_local[q].p = evaluate_polynomial(psi, eta, coeffs);

        //compute gradient:
        evaluate_grad(coeffs, psi, eta, h, &px, &py);
        
        //compute normal derivative of pressure at quadrature point:
        nx = Ctx->surface.q_normal[q].x;
        ny = Ctx->surface.q_normal[q].y;
        
        Ctx->surface.q_local[q].pr = (px/h)*nx + (py/h)*ny;
        
    }

    //gather Kirchhoff data to the 0th process and write it to a file:

    // as we are transfering struct... need a custum data type
    // 1- Here, create all the properties to call MPI_Type_create_struct
    MPI_Datatype kirchhoff_data_t;

    MPI_Aint displacements[3] = {offsetof(KirchhoffData, theta), offsetof(KirchhoffData, p), offsetof(KirchhoffData, pr)};
    int block_lengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

    // 2- Create the type, and commit it
    MPI_Type_create_struct(3, block_lengths, displacements, types, &kirchhoff_data_t);
    MPI_Type_commit(&kirchhoff_data_t);

    //size of vector:
    int local_size_per_proc[nproc];
    int local_size = Ctx->surface.Nk_local;

    ierr = MPI_Allgather(&local_size, 1, MPI_INT, local_size_per_proc, 1, MPI_INT, PETSC_COMM_WORLD); CHKERRQ(ierr);

     // displacement vector:
    int disp[nproc];
    int proc_begin = 0;
    for (k = 0; k < rank; ++k)
        proc_begin += local_size_per_proc[k];

    ierr = MPI_Allgather(&proc_begin, 1, MPI_INT, disp, 1, MPI_INT, MPI_COMM_WORLD); CHKERRQ(ierr);
    
    ierr = MPI_Gatherv(Ctx->surface.q_local, Ctx->surface.Nk_local, kirchhoff_data_t, Ctx->surface.q_global, local_size_per_proc, disp, kirchhoff_data_t, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);
    
    MPI_Type_free(&kirchhoff_data_t);

    return ierr;
}



//----------------------------------------------------------------------------
// Free memory for Kirchhoff data 
//----------------------------------------------------------------------------
PetscErrorCode FreeKirchhoffData(AppCtx* Ctx){

    PetscErrorCode ierr;
    int rank;

    free(Ctx->surface.q_local);
    free(Ctx->surface.q_point);
    free(Ctx->surface.q_normal);
    free(Ctx->surface.q_cell);
    
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr); 

    if (rank == 0)
        free(Ctx->surface.q_global);

    return ierr;
}