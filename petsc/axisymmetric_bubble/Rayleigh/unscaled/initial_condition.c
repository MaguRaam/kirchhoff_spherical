/*
 * initial_condition.c
 *      Author: sunder
 */
#include "hype.h"

//----------------------------------------------------------------------------
// Common Initial Conditions
//----------------------------------------------------------------------------

PetscReal TanhRadial(PetscReal, PetscReal, PetscReal, PetscReal, PetscReal, PetscReal, PetscReal, PetscReal);
void RayleighCollapse(PetscReal, PetscReal, PetscReal *);
void WallCollapse(PetscReal, PetscReal, PetscReal *);
void CavityCollapse(PetscReal, PetscReal, PetscReal *);

//----------------------------------------------------------------------------
// Initial condition function
//----------------------------------------------------------------------------

void InitialCondition(PetscReal x, PetscReal y, PetscReal *Q0)
{
    RayleighCollapse(x, y, Q0);
}

PetscReal TanhRadial(PetscReal x, PetscReal y, PetscReal x0, PetscReal y0, PetscReal R0, PetscReal eps, PetscReal in, PetscReal out)
{

    PetscReal R = PetscSqrtReal((x - x0) * (x - x0) + (y - y0) * (y - y0));

    return 0.5 * ((in + out) + (out - in) * PetscTanhReal((R - R0) / (eps)));
}

//----------------------------------------------------------------------------
// Rayleigh Collapse of a Bubble
// [r,z] \in [0,10R] x [-10R,10R] (R = 0.038)
// Final Time: 3Tc, Tc = 3.5e-4
// BC: L-T, R-T, B-R, T-R
// GAMMA_1 = 4.4; GAMMA_2 = 1.4; PI_1 = 6.0e8; PI_2 = 0
//----------------------------------------------------------------------------

void RayleighCollapse(PetscReal x, PetscReal y, PetscReal *Q0)
{

    PetscReal V0[nVar];

    PetscReal eps = 1.0e-5;
    PetscReal x0 = 0.0, y0 = 0.0; // Center of the bubble
    PetscReal smear = 1.0;
    PetscReal h = 10.0 * R0 / 250.0;

    PetscReal p_inf = 1.0e6;  // Far field Pressure
    PetscReal p_b = 1.0e5;    // Pressure inside the bubble

    PetscReal rho_inf = 1000.; // Far field density
    PetscReal rho_b = 1.0; // Bubble density

    V0[0] = TanhRadial(x, y, x0, y0, R0, smear * h, rho_b, rho_inf);
    V0[1] = 0.0;
    V0[2] = 0.0;
    V0[3] = TanhRadial(x, y, x0, y0, R0, smear * h, p_b, p_inf);
    V0[4] = TanhRadial(x, y, x0, y0, R0, smear * h, eps, 1.0 - eps);

    PDEPrim2Cons(V0, Q0);
}

//----------------------------------------------------------------------------
// Bubble Collapse Near Wall
// [x,y] \in [0,10] x [0.0,2.5]
// Final Time: 0.02
// BC: L-R, R-T, B-T, T-T
// g1 = 4.4; g2 = 1.4; p1 = 6000.0; g1 = 0.0
//----------------------------------------------------------------------------

void WallCollapse(PetscReal x, PetscReal y, PetscReal *Q0)
{

    PetscReal V0[nVar];

    PetscReal R0 = 1.0;
    PetscReal x0 = 0.0;
    const PetscReal y0 = 0.0;
    PetscReal R = PetscSqrtReal((x - x0) * (x - x0) + (y - y0) * (y - y0));

    // Inside the bubble
    PetscReal rhoi = 1.0e-2;
    PetscReal pi = 1.0;
    PetscReal alphai = 0.0;

    // Outside the bubble
    PetscReal rhoo = 1.0;
    PetscReal po = 100.0;
    PetscReal alphao = 1.0;

    PetscReal h = 5.0 / 500.0;
    PetscReal smear = 1.0;

    V0[0] = 0.5 * ((rhoi + rhoo) + (rhoo - rhoi) * PetscTanhReal((R - R0) / (smear * h)));
    V0[1] = 0.0;
    V0[2] = 0.0;
    V0[3] = 0.5 * ((pi + po) + (po - pi) * PetscTanhReal((R - R0) / (smear * h)));
    V0[4] = 0.5 * ((alphai + alphao) + (alphao - alphai) * PetscTanhReal((R - R0) / (smear * h)));

    PDEPrim2Cons(V0, Q0);
}

//----------------------------------------------------------------------------
// Mach 1.72 Shock in water hitting air cavity
// [x,y] \in [0,12] x [0,6]
// Final Time: 0.045
// BC: L-T, R-T, B-T, T-T
// g1 = 4.4; g2 = 1.4; p1 = 6000.0; g1 = 0.0
//----------------------------------------------------------------------------

void CavityCollapse(PetscReal x, PetscReal y, PetscReal *Q0)
{

    PetscReal V0[nVar];

    PetscReal x0 = 6.0;
    PetscReal y0 = 0.0;
    PetscReal R0 = 3.0;
    PetscReal eps = 1.0e-5;
    PetscReal h = 12.0 / 1024.0;
    PetscReal smear = 1.0;

    V0[0] = TanhRadial(x, y, x0, y0, R0, smear * h, 0.0012, 1.0);
    V0[1] = 0.0;
    V0[2] = 0.0;
    V0[3] = 1.0;
    V0[4] = TanhRadial(x, y, x0, y0, R0, smear * h, eps, 1.0 - eps);

    // (Water) Shock

    if (x > 11.4)
    {

        V0[0] = 1.325;
        V0[1] = -68.525;
        V0[2] = 0.0;
        V0[3] = 19153.0;
        V0[4] = 1.0 - eps;
    }

    PDEPrim2Cons(V0, Q0);
}
