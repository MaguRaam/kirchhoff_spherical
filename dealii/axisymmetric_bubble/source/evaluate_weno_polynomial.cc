#include "../include/Weno432.h"


// Evaluate the WENO polynomial at a point 


double evaluate_weno_polynomial(const Vector<double>& coeffs, const Vector<double>& consts, const Point<2>& P, double h) {
    
	double dx = (P(0)-consts(0))/h; double dy = (P(1)-consts(1))/h; 
    
    // Use Honers Algorith 
    
    double retval = coeffs(0) -  coeffs(3)*consts(2) - coeffs(4)*consts(3) - coeffs(5)*consts(4) - 
			        coeffs(6)*consts(5) - coeffs(7)*consts(6) - coeffs(8)*consts(7) - coeffs(9)*consts(8) +
                    dx*(coeffs(1) + dy*coeffs(5) + dx*(coeffs(3) + dy*coeffs(8) + dx*coeffs(6)) ) + 
                    dy*(coeffs(2) + dy*(coeffs(4) + dy*coeffs(7) + dx*coeffs(9) ) ); 

    return retval; 
}

Point<2, double> evaluate_gradient(const Vector<double>& coeffs, const Vector<double>& consts, Point<2> P, double h) {
//void evaluate_weno_polynomial_grad(const Vector<double>& coeffs, const Vector<double>& consts,  const Point<2>& P, double* grad, const double& h) {

	Point<2, double> grad;
    double x = P(0); double y = P(1);
	
	double dx = (x-consts(0))/h; double dy = (y-consts(1))/h;
	double h_inv = 1.0/h;

	// Gradient in x direction
    grad[0] = h_inv * (coeffs(1) + 2.0*dx*coeffs(3) + dy*coeffs(5) + 3.0*dx*dx*coeffs(6) 
					+ 2.0*dx*dy*coeffs(8) + dy*dy*coeffs(9) );

	// Gradient in y direction
    grad[1] = h_inv * (coeffs(2) + 2.0*dy*coeffs(4) + dx*coeffs(5) + 3.0*dy*dy*coeffs(7) 
					+ dx*dx*coeffs(8) + 2.0*dx*dy*coeffs(9) );
					
    return grad; 					
}
