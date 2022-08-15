#include "../include/Weno432.h"

// Set the initial condtions  


double TanhRadial(double x,double y,double x0, double y0, double R0, double eps, double in, double out) {

    double R = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0));

    return 0.5*((in+out) + (out-in)*tanh((R-R0)/(eps)));
}


Vector<double> initial_condition(Point <2> p, double h) {


	Vector<double> W(5);

	
  // bubble radius:
  double R0 = 0.25;

	double eps = 1.0e-5;
  double x0 = 0.0, y0 = 0.0; // Center of the bubble
  double smear = 1.0;
  
  //scaled by 1e5
  double p_inf = 10;  			// Far field Pressure
  double p_b = 1;     			// Pressure inside the bubble

	//scaled by 1e3
	double rho_inf = 1.0;			// Far field density
	double rho_b = 0.001;			// Bubble density
  
  double x = p(0);
  double y = p(1);

  W[0] = TanhRadial(x,y,x0,y0,R0,smear*h,rho_b,rho_inf);
  W[1] = 0.0;
  W[2] = 0.0;
  W[3] = TanhRadial(x,y,x0,y0,R0,smear*h,p_b,p_inf);
  W[4] = TanhRadial(x,y,x0,y0,R0,smear*h,eps,1.0-eps);

	return W;
}


