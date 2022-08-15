#include "kirchhoff.h"


//Monopole solution:
class SphericalWave{

public:
    
    //constructor:
    SphericalWave(double c0, double f0) : c0(c0), f0(f0), t0(4.0/f0){}

    //pressure p
    double p(double r, double t){

		return source(t - r/c0)/(4.0*M_PI*r);
    }

    //pressure derivative dp/dr:
    double dpdr(double r, double t){
        
		return source_derivative(t - r/c0)/(-c0*4.0*M_PI*r) + source(t - r/c0)*(-0.25/(M_PI*r*r));
    }

    //pressure derivative dp/dt:
    double dpdt(double r, double t){
        
		return source_derivative(t - r/c0)/(4.0*M_PI*r);
    }

private:

	double source(double t){return -2.0*(t - t0)*f0*f0*exp( -1.0*f0*f0*(t - t0)*(t - t0));}
	double source_derivative(double t){return -1.0*f0*f0*f0*f0*(-2.0*t + 2.0*t0)*(2*t - 2*t0)*exp(-1.0*f0*f0*(t - t0)*(t - t0)) - 2.0*f0*f0*exp(-1.0*f0*f0*(t - t0)*(t - t0));}
	
	const double c0;
    const double f0;
	const double t0;
};



int main()
{

    // wave speed:
    const double c = 250.0;

    // radius of sphere:
    const double R = 1.0;

    // observer point:
    const Point xo{3.0, 3.0, 3.0};

    // observer time:
    double to = 0.0;
    const double tf = 0.1;
    const double dt = 0.001;

    // no of cells along theta and phi direction:
    const int n_theta = 100;
    const int n_phi = 2 * n_theta;

    // create theta and phi grid points:
    auto theta = create_grid(0, M_PI, n_theta);
    auto phi = create_grid(0, 2. * M_PI, n_phi);

    // compute quadrature points on sphere:
    auto xq = create_sphere(R, theta, phi);

	  // write sphere in vtk:
	  write_sphere(xq);

    // load pressure and time data:
    auto [t, p, pr] = load_pressure(n_theta);

    // compute distance and angle between quadrature point/normal and observer point:
    auto [rq, cosq] = compute_distance_angle(xq, xo);

    // write Kirchhoff data:
    std::ofstream observer("pressure.dat", std::ios::out);
    observer.flags(std::ios::dec | std::ios::scientific);
    observer.precision(16);
    
    //exact solution
  	SphericalWave wave(c, 100.);

    // loop over observer time and compute Kirchhoff integral
    while (to < tf)
    {
        std::cout << "observer time to = " << to << std::endl;

        auto tauq = compute_emission_time(xq, xo, to, c);

        observer << to << "\t" << compute_kirchhoff_integral(rq, cosq, tauq, theta, t, p, pr, R, c) << "\t"<< wave.p(norm(xo), to) << std::endl;

        to += dt;
    }

    return 0;
}
