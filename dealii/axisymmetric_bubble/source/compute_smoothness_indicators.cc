#include "../include/Weno432.h"

// Find smoothness indicators 

// Second order polynomial 

double compute_second_order_smoothness_indicator(const Vector<double>& coeffs, const Vector<double>& IS, double h) {

    double u_x   = coeffs(0); double u_y   = coeffs(1); 
    
    double const_term = IS(0)*(1./h)*(u_x*u_x + u_y*u_y);

    return (const_term);
	
}

double compute_second_order_smoothness_indicator(const Vector<double>& coeffs) {

    double u_x   = coeffs(0); double u_y   = coeffs(1); 
    
    double const_term = u_x*u_x + u_y*u_y;    
    return (const_term);
	
}

// Third Order polynomial 

double compute_third_order_smoothness_indicator(const Vector<double>& coeffs, const Vector<double>& IS, double h) {
    
    double u_x   = coeffs(0); double u_y   = coeffs(1); double u_xx  = coeffs(2); double u_yy  = coeffs(3);
    double u_xy  = coeffs(4); 

	double h_inv = 1.0/h;
	
    double const_term = IS(0)*(h_inv)*(4*u_xx*u_xx + u_xy*u_xy + 4*u_yy*u_yy + u_x*u_x + u_y*u_y);

    double x_term = IS(1)*(h_inv)*(4*u_x*u_xx + 2*u_xy*u_y);

    double y_term = IS(2)*(h_inv)*(2*u_x*u_xy + 4*u_y*u_yy); 

    double x2_term = IS(3)*(h_inv)*(4*u_xx*u_xx + u_xy*u_xy);

    double y2_term = IS(4)*(h_inv)*(u_xy*u_xy + 4*u_yy*u_yy); 

    double xy_term = IS(5)*(h_inv)*(4*u_xx*u_xy + 4*u_xy*u_yy);
    
    return (const_term + x_term + y_term + x2_term + y2_term + xy_term); 

	
	
	/*
	
	double exact = h*(6*u_x*u_x + 26*u_xx*u_xx + 7*u_xy*u_xy + 6*u_y*u_y + 26*u_yy*u_yy)/6.;

	return exact; 

	*/ 
}

double compute_third_order_smoothness_indicator(const Vector<double>& coeffs) {
    
    double u_x   = coeffs(0); double u_y   = coeffs(1); double u_xx  = coeffs(2); double u_yy  = coeffs(3);
    double u_xy  = coeffs(4); 

    double const_term = u_xx*u_xx + u_xy*u_xy + u_yy*u_yy + u_x*u_x + u_y*u_y;	

    return const_term;
 
}

// Fourth order polynomial 

double compute_fourth_order_smoothness_indicator(const Vector<double>& coeffs, const Vector<double>& IS, double h) {
    
    double u_x   = coeffs(0); double u_y   = coeffs(1); double u_xx  = coeffs(2); double u_yy  = coeffs(3);
    double u_xy  = coeffs(4); double u_xxx = coeffs(5); double u_yyy = coeffs(6); double u_xxy = coeffs(7);
    double u_xyy = coeffs(8);
	
	double h_inv = 1.0/h;
	
	double const_term = IS(0)*(h_inv)*(36*u_xxx*u_xxx + 4*u_xxy*u_xxy + 4*u_xyy*u_xyy + 36*u_yyy*u_yyy + 
                               4*u_xx*u_xx + u_xy*u_xy + 4*u_yy*u_yy + u_x*u_x + u_y*u_y);

    double x_term = IS(1)*(h_inv)*(24*u_xx*u_xxx + 4*u_xxy*u_xy + 8*u_xyy*u_yy + 4*u_x*u_xx + 2*u_xy*u_y);

    double y_term = IS(2)*(h_inv)*(8*u_xx*u_xxy + 4*u_xy*u_xyy + 24*u_yy*u_yyy + 2*u_x*u_xy + 4*u_y*u_yy); 

    double x2_term = IS(3)*(h_inv)*(36*u_xxx*u_xxx + 4*u_xxy*u_xxy + 4*u_xyy*u_xyy + 6*u_x*u_xxx + 
                             4*u_xx*u_xx + 2*u_xxy*u_y + u_xy*u_xy);

    double y2_term = IS(4)*(h_inv)*(4*u_xxy*u_xxy + 4*u_xyy*u_xyy + 36*u_yyy*u_yyy + 
                            2*u_x*u_xyy + u_xy*u_xy + 6*u_y*u_yyy + 4*u_yy*u_yy); 

    double xy_term = IS(5)*(h_inv)*(24*u_xxx*u_xxy + 8*u_xxy*u_xyy + 24*u_xyy*u_yyy + 4*u_x*u_xxy + 
                             4*u_xx*u_xy + 4*u_xy*u_yy + 4*u_xyy*u_y);

    double x3_term = IS(6)*(h_inv)*(12*u_xx*u_xxx + 2*u_xxy*u_xy); 

    double y3_term = IS(7)*(h_inv)*(2*u_xy*u_xyy + 12*u_yy*u_yyy); 

    double x2y_term = IS(8)*(h_inv)*(8*u_xx*u_xxy + 6*u_xxx*u_xy + 4*u_xxy*u_yy + 4*u_xy*u_xyy);

    double xy2_term = IS(9)*(h_inv)*(4*u_xx*u_xyy + 4*u_xxy*u_xy + 6*u_xy*u_yyy + 8*u_xyy*u_yy); 

    double x4_term = IS(10)*(h_inv)*(9*u_xxx*u_xxx + u_xxy*u_xxy);

    double y4_term = IS(11)*(h_inv)*(u_xyy*u_xyy + 9*u_yyy*u_yyy); 

    double x2y2_term = IS(12)*(h_inv)*(6*u_xxx*u_xyy + 4*u_xxy*u_xxy + 6*u_xxy*u_yyy + 4*u_xyy*u_xyy); 

    double x3y_term = IS(13)*(h_inv)*(12*u_xxx*u_xxy + 4*u_xxy*u_xyy); 

    double xy3_term = IS(14)*(h_inv)*(4*u_xxy*u_xyy + 12*u_xyy*u_yyy);

    
    return (const_term + x_term + y_term + x2_term + y2_term + xy_term + x3_term + 
           y3_term + x2y_term + xy2_term + x4_term + y4_term + x2y2_term + x3y_term + xy3_term); 
	
} 

double compute_fourth_order_smoothness_indicator(const Vector<double>& coeffs) {
    
    double u_x   = coeffs(0); double u_y   = coeffs(1); double u_xx  = coeffs(2); double u_yy  = coeffs(3);
    double u_xy  = coeffs(4); double u_xxx = coeffs(5); double u_yyy = coeffs(6); double u_xxy = coeffs(7);
    double u_xyy = coeffs(8);

	double const_term = u_xxx*u_xxx + u_xxy*u_xxy + u_xyy*u_xyy + u_yyy*u_yyy + 
                               u_xx*u_xx + u_xy*u_xy + u_yy*u_yy + u_x*u_x + u_y*u_y;

    return const_term;

} 
