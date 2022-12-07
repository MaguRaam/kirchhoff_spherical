/*
 * post_process.c
 *      Author: sunder
 */ 

#include "hype.h"  



//----------------------------------------------------------------------------
// Find total volume fraction in the domain
//----------------------------------------------------------------------------

PetscErrorCode PP(Vec U, DM da, PetscReal h, PetscReal* vol, PetscReal* p_max, PetscReal* tau_xx, PetscReal *tau_yy, PetscReal* tau_zz, PetscReal* tau_xy) {

    PetscErrorCode ierr=0;

    Field  **u;
    PetscInt i, j, xs, ys, xm, ym, c;
    PetscReal vol_local = 0.0;
    PetscReal vol_global;
    PetscReal h2 = h*h;
    PetscReal pmax = 0.0, xxmax = 0.0, yymax = 0.0, zzmax = 0.0, xymax = 0.0;
    PetscReal y, xx, yy, xy, zz;
    PetscReal phi, u_x, v_x, u_y, v_y, prs, mu, v, div_v;
    PetscReal Q[nVar], V[nVar];

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, U, &u);CHKERRQ(ierr);


    for (j = ys; j < ys+ym; ++j) {
        for (i = xs; i < xs+xm; ++i) {

            phi = u[j][i].comp[4];
            vol_local += h2*phi;

            if (i == 0) {

                y = ((PetscReal)(j)+0.5)*h;

                for (c = 0; c < nVar; ++c)
                    Q[c] = u[j][i].comp[c];


                PDECons2Prim(Q,V);

                prs = V[3];

                mu = phi*mu1 + (1.0-phi)*mu2;
                v = u[j][i].comp[2]/u[j][i].comp[0];


                u_x = (u[j][i+1].comp[1]/u[j][i+1].comp[0] - u[j][i].comp[1]/u[j][i].comp[0])/h;
                v_x = (u[j][i+1].comp[2]/u[j][i+1].comp[0] - u[j][i].comp[2]/u[j][i].comp[0])/h;

                if (j == ys) {

                    u_y = (u[j+1][i].comp[1]/u[j+1][i].comp[0] - u[j][i].comp[1]/u[j][i].comp[0])/h;
                    v_y = (u[j+1][i].comp[2]/u[j+1][i].comp[0] - u[j][i].comp[2]/u[j][i].comp[0])/h;

                }

                else if (j == ys+ym-1) {

                    u_y = (u[j][i].comp[1]/u[j][i].comp[0] - u[j-1][i].comp[1]/u[j-1][i].comp[0])/h;
                    v_y = (u[j][i].comp[2]/u[j][i].comp[0] - u[j-1][i].comp[2]/u[j-1][i].comp[0])/h;

                }

                else {

                    u_y = 0.5*(u[j+1][i].comp[1]/u[j+1][i].comp[0] - u[j-1][i].comp[1]/u[j-1][i].comp[0])/h;
                    v_y = 0.5*(u[j+1][i].comp[2]/u[j+1][i].comp[0] - u[j-1][i].comp[2]/u[j-1][i].comp[0])/h;

                }


                div_v = u_x + v_y + v/y;


                xx = mu*(2.0*u_x - (2./3.)*div_v);
                xy = mu*(u_y + v_x);
                yy = mu*(2.0*v_y - (2./3.)*div_v);
                zz = mu*(2.0*v/y - (2./3.)*div_v);

                if (prs > pmax)
                    pmax = prs;

                if (PetscAbsReal(xx) > xxmax)
                    xxmax = xx;


                if (PetscAbsReal(yy) > yymax)
                    yymax = yy;


                if (PetscAbsReal(zz) > zzmax)
                    zzmax = zz;


                if (PetscAbsReal(xy) > xymax)
                    xymax = xy;
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(da, U, &u);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&vol_local, &vol_global,1,MPIU_REAL,MPIU_SUM, PETSC_COMM_WORLD);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&pmax, p_max,1,MPIU_REAL,MPIU_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&xxmax, tau_xx, 1,MPIU_REAL,MPIU_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&yymax, tau_yy, 1,MPIU_REAL,MPIU_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&zzmax, tau_zz, 1,MPIU_REAL,MPIU_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);
    ierr = MPI_Allreduce(&xymax, tau_xy, 1,MPIU_REAL,MPIU_MAX, PETSC_COMM_WORLD);CHKERRQ(ierr);

    *vol = vol_global;

    return ierr;
}

//----------------------------------------------------------------------------
// Monitor function for additional processing in the intermediate time steps 
//----------------------------------------------------------------------------

PetscErrorCode MonitorFunction (TS ts,PetscInt step, PetscReal time, Vec U, void *ctx) {

    PetscErrorCode ierr;
    AppCtx *Ctx = (AppCtx*)ctx;
    DM da;
    PetscMPIInt MyPID;
    PetscReal vol, prs, tau_xx, tau_yy, tau_zz, tau_xy;

    ierr = TSGetDM(ts,&da);CHKERRQ(ierr);
    // Get rank of the processor
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&MyPID);CHKERRQ(ierr);
    // Set the time step based on CFL condition
    ierr = TSSetTimeStep(ts, Ctx->dt);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%d t = %.5e\n", step, time);CHKERRQ(ierr);

    // Plot the solution at the required time interval

    if (Ctx->WriteInterval != 0) {

        if(step%Ctx->WriteInterval == 0) {

            ierr = ComputePrimitiveVariables(U, Ctx->W, da);CHKERRQ(ierr);
            char filename[30];
            sprintf(filename, "plot/sol-%08d.vts", step); // 8 is the padding level, increase it for longer simulations
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data in vts format to %s at t = %f, step = %d\n", filename, time, step);CHKERRQ(ierr);
            PetscViewer viewer;
            ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);
            ierr = DMView(da, viewer);
            VecView(Ctx->W, viewer);
            ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

        }
    }

    if (Ctx->RestartInterval != 0) {


        PetscViewer viewer_binary;

        if(step%Ctx->RestartInterval == 0) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data in binary to restart1.bin at t = %f\n", time);CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"restart1.bin",FILE_MODE_WRITE, &viewer_binary);CHKERRQ(ierr);
            ierr = VecView(U,viewer_binary);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer_binary); CHKERRQ(ierr);
        }

        if((step+10)%(Ctx->RestartInterval) == 0) {
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing data in binary to restart2.bin at t = %f\n", time);CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"restart2.bin",FILE_MODE_WRITE, &viewer_binary);CHKERRQ(ierr);
            ierr = VecView(U,viewer_binary);CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer_binary);CHKERRQ(ierr);
        }
    }


    ierr= PP(U, da, Ctx->h, &vol, &prs, &tau_xx, &tau_yy, &tau_zz, &tau_xy);CHKERRQ(ierr);
    if (MyPID == 0) {
        FILE* fp;
        fp = fopen("WallData.csv", "a");
        fprintf(fp, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n", time, vol, prs, tau_xx,tau_yy,tau_zz,tau_xy);
        fclose(fp);
    }

    //write Kirchhoff data:

    if (step%Ctx->KirchhoffWriteInterval == 0){

        //compute p', p'r on Kirchhoff surface:			
        ierr =  WriteKirchhoffData(U, da, Ctx, time); CHKERRQ(ierr);


        if (MyPID == 0) 
        {
    
            FILE *t_file = fopen("kirchhoff/data/time.dat","a+");
            FILE *p_file = fopen("kirchhoff/data/pressure.dat","a+");

            //write time.dat
            ierr = PetscFPrintf(PETSC_COMM_SELF, t_file, "%.16e\n", time); CHKERRQ(ierr);

            //loop over quadrature points and write pressure data:
            ierr = PetscFPrintf(PETSC_COMM_SELF, p_file, "t = %.16e\n", time); CHKERRQ(ierr);

            for (PetscInt q = 0; q < Ctx->surface.Nk; ++q){
            	ierr = PetscFPrintf(PETSC_COMM_SELF, p_file, "%.16e\t%.16e\t%.16e\n", Ctx->surface.q_global[q].theta, Ctx->surface.q_global[q].p, Ctx->surface.q_global[q].pr); CHKERRQ(ierr); 

            }
            fclose(p_file); 
            fclose(t_file);

        }
    
    }       

    return ierr;
}

//----------------------------------------------------------------------------
// Find L2 and L_inf errors for periodic test cases with square domain and 
// where final solution conicides with the initial condition
//----------------------------------------------------------------------------

PetscErrorCode ErrorNorms(Vec U, DM da, AppCtx Ctx, PetscReal* l2, PetscReal* linf) {

    PetscErrorCode ierr;

    DM          coordDA;
    Vec         coordinates;
    DMDACoor2d  **coords;
    Field   **u;
    Field   **u_exact;
    PetscInt    xs, ys, xm, ym, i, j, c, l , m;
    PetscReal integral[nVar]; 
    PetscReal xc, yc, xGP, yGP;
    PetscReal h = Ctx.h;
    PetscReal Q0[nVar];
    Vec U_exact;
    
    ierr = VecDuplicate(U,&U_exact);CHKERRQ(ierr);
    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);CHKERRQ(ierr);
    ierr = DMGetCoordinateDM(da, &coordDA);CHKERRQ(ierr);
    ierr = DMGetCoordinates(da, &coordinates);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(coordDA, coordinates, &coords);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, U_exact, &u_exact);CHKERRQ(ierr);

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
            
            for (c = 0; c < nVar; ++c) {
                if (c == 0) {
                    u_exact[j][i].comp[c] = integral[c]; 
                }
                
                else {
                    u[j][i].comp[c] = 0.0;
                    u_exact[j][i].comp[c] = 0.0; 
                }
            }
        }
    }

    ierr = DMDAVecRestoreArray(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, U_exact, &u_exact);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(coordDA, coordinates, &coords);CHKERRQ(ierr);

    PetscReal nrm_inf, nrm2; 

    ierr = VecAXPY(U_exact, -1.0, U);CHKERRQ(ierr);
    ierr = VecNorm(U_exact, NORM_INFINITY, &nrm_inf);CHKERRQ(ierr);
    ierr = VecNorm(U_exact, NORM_2, &nrm2);CHKERRQ(ierr);

    nrm2 = nrm2/((PetscReal)Ctx.N_x);

    *l2 = nrm2; 
    *linf = nrm_inf; 

    ierr = VecDestroy(&U_exact);CHKERRQ(ierr);

    return ierr; 
}

//----------------------------------------------------------------------------
// Find the cell averages of primitive variables in each cell 
//----------------------------------------------------------------------------

PetscErrorCode ComputePrimitiveVariables(Vec U, Vec W, DM da) {
    
    PetscErrorCode ierr;           // For catching PETSc errors 
    PetscInt c,i,j,xs,ys,xm,ym;    // Corners of the grid on the given solution 
    Field  **u;                    // Local array of the conserved variables 
    Field  **w;                    // Local array of the primitive variables 
    PetscReal Q[nVar], V[nVar];  

    // Read the local solution to the array u  

    ierr = DMDAVecGetArrayRead(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, W, &w);CHKERRQ(ierr);

    ierr = DMDAGetCorners(da, &xs, &ys, NULL, &xm, &ym, NULL);


    for (j=ys; j<ys+ym; ++j) {
        for (i=xs; i<xs+xm; ++i) {
        
            // Select stencils for each component of the solution 
    
            for (c = 0 ; c < nVar; ++c)
                Q[c] = u[j][i].comp[c];
    
            PDECons2Prim(Q, V);
            
            for (c = 0 ; c < nVar; ++c)
                w[j][i].comp[c] = V[c];
        }
    }  

    ierr = DMDAVecRestoreArrayRead(da, U, &u);CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da,W,&w);CHKERRQ(ierr);

    return ierr; 
}

