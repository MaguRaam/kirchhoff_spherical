#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../include/trueSolution.h"

Point<2> grid_transform (Point<2>& in)
{

	double x = in(0), x_ = 0.0, x_max = 0.0, xc = 0.0, lengthx = 10.0, zeta0 = 0.55;
	double A0 = 8.0;
//	double A0 = 10.0;
	double y_ = 0.0, y_max = 0.0, yc = 0.0, lengthy = 2.5 ;
	double y = std::fabs(in(1));


	double zeta_x = (x-xc)/double(lengthx) ;
	double zeta_y = (y-yc)/double(lengthy) ;
	x_max = zeta0 + sinh(A0*(1.0 - zeta0))/A0;
	y_max = zeta0 + sinh(A0*(1.0 - zeta0))/A0;

	if (zeta_x > zeta0)
		x_ = zeta0 + sinh(A0*(zeta_x - zeta0))/A0;
	else 
		x_ = zeta_x;

	if (zeta_y > zeta0)
		y_ = zeta0 + sinh(A0*(zeta_y - zeta0))/A0;
	else 
		y_ = zeta_y;

	x_ = x_*lengthx/x_max;
	y_ = y_*lengthy/y_max;

	if(in(1) < 0 ) y_ = -y_;
		
	Point<2> out(x_, y_);

//	if (y == -1 || x == 0) std::cout<<"input: "<<in<<"\tout: "<<out<<"\tx_max: "<<x_max<<std::endl;
	return out;
}
/*
Point<2>
SineTransform::push_forward(const Point<2> &chart_point) const
{
	double x = chart_point(0), y = chart_point(1);
	double perturb_y = 0.0, perturb_y0 = 0.0, perturb_y1 = 0.0, perturb_y2 = 0.0, perturb_y3 = 0.0, perturb_y4 = 0.0, 
			x0 = 4.5, x1 = 4.75, x2 = 5.0, x3 = 5.25, x4 = 5.5, sigma = 0.075;

//	if(x >= 2.0 && x <= 4.0)
//		perturb_y = 0.08*exp(-10.0*(x - x1)*(x - x1));


	perturb_y0 = 0.08*exp(-((x - x0)/sigma)*((x - x0)/sigma));
	perturb_y1 = 0.08*exp(-((x - x1)/sigma)*((x - x1)/sigma));
	perturb_y2 = 0.08*exp(-((x - x2)/sigma)*((x - x2)/sigma));
	perturb_y3 = 0.08*exp(-((x - x3)/sigma)*((x - x3)/sigma));
	perturb_y4 = 0.08*exp(-((x - x4)/sigma)*((x - x4)/sigma));

	perturb_y = perturb_y0 + perturb_y1 + perturb_y2 + perturb_y3 + perturb_y4;

//	perturb_y = 0.08*exp(-((x - x1)/sigma)*((x - x1)/sigma));

	Point<2> out(chart_point(0), chart_point(1) + perturb_y);
	return out;
}

Point<2>
SineTransform::pull_back(const Point<2> &space_point) const
{

	double x = space_point(0), y = space_point(1);
	double perturb_y = 0.0, perturb_y0 = 0.0, perturb_y1 = 0.0, perturb_y2 = 0.0, perturb_y3 = 0.0, perturb_y4 = 0.0, 
			x0 = 4.5, x1 = 4.75, x2 = 5.0, x3 = 5.25, x4 = 5.5, sigma = 0.075;

//	if(x >= 2.0 && x <= 4.0)
//		perturb_y = 0.08*exp(-10.0*(x - x1)*(x - x1));


//	perturb_y = 0.08*exp(-((x - x1)/sigma)*((x - x1)/sigma));

	perturb_y0 = 0.08*exp(-((x - x0)/sigma)*((x - x0)/sigma));
	perturb_y1 = 0.08*exp(-((x - x1)/sigma)*((x - x1)/sigma));
	perturb_y2 = 0.08*exp(-((x - x2)/sigma)*((x - x2)/sigma));
	perturb_y3 = 0.08*exp(-((x - x3)/sigma)*((x - x3)/sigma));
	perturb_y4 = 0.08*exp(-((x - x4)/sigma)*((x - x4)/sigma));

	perturb_y = perturb_y0 + perturb_y1 + perturb_y2 + perturb_y3 + perturb_y4;

	Point<2> out(space_point(0), space_point(1) - perturb_y);
	return out;

}


std::unique_ptr<Manifold<2,2> >
SineTransform::clone() const
{
  return std_cxx14::make_unique<SineTransform>();
}
*/
