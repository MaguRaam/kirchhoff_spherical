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
// Compute no of quadpoints (on Kirchhoff arc) in the current process: 
//----------------------------------------------------------------------------
PetscErrorCode ComputeNofQuadPoints(DM da, AppCtx* Ctx){

    PetscErrorCode  ierr;
    PetscInt        xs, ys, xm, ym, i, j, k;
    PetscReal       theta, xq, yq;
    PetscReal       Rk = Ctx->surface.Rk;
    PetscInt        Nk = Ctx->surface.Nk;

    // compute dtheta:
    Ctx->surface.dtheta = PETSC_PI / (PetscReal)Nk;
    PetscReal dtheta = Ctx->surface.dtheta;

    // loop over quad points on Kirchhoff arc:
    Ctx->surface.Nk_local = 0;

    // local domain boundaries:
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    for (k = 0; k < Nk; ++k)
    {
        //compute theta, x, y at quadrature point:
        theta = (k + 0.5) * dtheta;
        xq = Rk * PetscCosReal(theta); // z
        yq = Rk * PetscSinReal(theta); // r

        //get cell index of the qudrature point:
        i = GetCellIndex(xq, Ctx->x_min, Ctx->h);
        j = GetCellIndex(yq, Ctx->y_min, Ctx->h);

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
    PetscReal theta, xq, yq, distanceq;
    PetscReal dtheta = Ctx->surface.dtheta;
    int rank;

    //allocate memory for local arary:
    Ctx->surface.q_local = (KirchhoffData*)malloc(Nk_local*sizeof(KirchhoffData));
    Ctx->surface.q_point = (Point*)malloc(Nk_local*sizeof(Point));
    Ctx->surface.q_normal = (Point*)malloc(Nk_local*sizeof(Point));
    Ctx->surface.q_cell  = (CellIndex*)malloc(Nk_local*sizeof(CellIndex));     

    //allocate memory in rank 0 for global data:
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr); 
    
    if (rank == 0)
        Ctx->surface.q_global = (KirchhoffData*)malloc(Nk*sizeof(KirchhoffData));

     // local domain boundaries:
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);

    //loop over quadarture points in the current process:
    int q = 0;
    for (k = 0; k < Nk; ++k){

        //compute quadrature point:
        theta = (k + 0.5) * dtheta;
        xq = Rk * PetscCosReal(theta); // z
        yq = Rk * PetscSinReal(theta); // r

        //get cell index of the qudrature point:
        i = GetCellIndex(xq, Ctx->x_min, Ctx->h);
        j = GetCellIndex(yq, Ctx->y_min, Ctx->h);

        //check if the point is in the current process:
        if (i >= xs &&  i < xs+xm && j >= ys &&  j < ys+ym){

            //store quadrature point:
            Ctx->surface.q_point[q].x = xq;
            Ctx->surface.q_point[q].y = yq;

            //compute distance from center for quadrature point:
            distanceq = PetscSqrtReal(xq*xq + yq*yq);

            //store normal vector:
            Ctx->surface.q_normal[q].x = xq/distanceq;
            Ctx->surface.q_normal[q].y = yq/distanceq;

            //store cell index of quadrature point:
            Ctx->surface.q_cell[q].i = i;
            Ctx->surface.q_cell[q].j = j;

            //store theta:
            Ctx->surface.q_local[q].theta = theta;

            q++;
        }

    }

    return ierr;
}

//----------------------------------------------------------------------------------------------------
// Get perturbed pressure p' = p - p0 from conserved state variable
//----------------------------------------------------------------------------------------------------
PetscReal GetPerturbedPressure(const Field* q){

    PetscReal Q[nVar], V[nVar];

    for (int c = 0 ; c < nVar; ++c) Q[c] = q->comp[c];
    PDECons2Prim(Q, V);

    return V[3] - p0;
}

//----------------------------------------------------------------------------------------------------
// Reconstruct polynomial on a given cell
//----------------------------------------------------------------------------------------------------
void ReconstructOnGivenCell(Field** u, PetscInt i, PetscInt j, PetscReal coeffs[], PetscReal pcoeffs[]){

    PetscReal p_x_loc[s_width], p_y_loc[s_width], p_xy_loc[s_width];

    //get perturbed pressure values from the stencil:
    p_x_loc[0] = GetPerturbedPressure(&u[j][i-2]); 
    p_x_loc[1] = GetPerturbedPressure(&u[j][i-1]); 
    p_x_loc[2] = GetPerturbedPressure(&u[j][i]); 
    p_x_loc[3] = GetPerturbedPressure(&u[j][i+1]); 
    p_x_loc[4] = GetPerturbedPressure(&u[j][i+2]);

    p_y_loc[0] = GetPerturbedPressure(&u[j-2][i]); 
    p_y_loc[1] = GetPerturbedPressure(&u[j-1][i]); 
    p_y_loc[2] = GetPerturbedPressure(&u[j][i]); 
    p_y_loc[3] = GetPerturbedPressure(&u[j+1][i]); 
    p_y_loc[4] = GetPerturbedPressure(&u[j+2][i]);

    p_xy_loc[0] = GetPerturbedPressure(&u[j][i]);
    p_xy_loc[1] = GetPerturbedPressure(&u[j+1][i+1]);
    p_xy_loc[2] = GetPerturbedPressure(&u[j-1][i+1]);
    p_xy_loc[3] = GetPerturbedPressure(&u[j+1][i-1]);
    p_xy_loc[4] = GetPerturbedPressure(&u[j-1][i-1]);

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
    PetscReal xc, yc, psi, eta, px, py, nx, ny, h = Ctx->h;
    Field **u;

    // get rank and no of process:
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&nproc);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr); 

    // Scatter global->local to have access to the required ghost values 
    ierr = DMGlobalToLocalBegin(da, U, INSERT_VALUES, Ctx->localU);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(da, U, INSERT_VALUES,   Ctx->localU);CHKERRQ(ierr);

    // Read the local solution to the array u  
    ierr = DMDAVecGetArrayRead(da, Ctx->localU, &u); CHKERRQ(ierr); 

    for (q = 0; q < Ctx->surface.Nk_local; ++q){    
        
        //compute cell index containing quad point:
        i = Ctx->surface.q_cell[q].i;
        j = Ctx->surface.q_cell[q].j;

        //reconstruct polynomial on a cell containing quadrature point:
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

    //! as we are transfering struct... need a custum data type
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
