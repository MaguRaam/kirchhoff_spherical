#ifndef RIEMANN_SOLVERS_H_
#define RIEMANN_SOLVERS_H_

#include "Headers.h"

// Function Declarations

Vector<double> primitive_to_conserved(Vector<double>);
Vector<double> conserved_to_primitive(Vector<double>);
Vector<double> compute_flux_from_primitive_variable(Vector<double>, double, double);
Vector<double> local_Lax_Friedrichs_riemann_solver(Vector<double>, Vector<double>, double, double, Point<2>, bool);
Vector<double> HLLC_riemann_solver(Vector<double>, Vector<double>, double, double, Point<2>, bool);
Vector<double> rotated_HLLC_riemann_solver(Vector<double>, Vector<double>, double, double, Point<2>, bool);



#endif /* RIEMANN_SOLVERS_H_ */
