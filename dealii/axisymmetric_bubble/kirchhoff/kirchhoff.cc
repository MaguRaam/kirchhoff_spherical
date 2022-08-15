#include "kirchhoff.h"


int main()
{

    // wave speed:
    const double c = 162.61611236282832;

    // radius of sphere:
    const double R = 0.25;

    // observer point:
    const Point xo{0.0, 29.0*R, 0.0};

    // observer time:
    double to = 0.1;
    const double tf = 0.15;
    const double dt = 0.001;

    // no of cells along theta and phi direction:
    const int n_theta = 5000;
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
    std::ofstream observer("p_kirchhoff.dat", std::ios::out);
    observer.flags(std::ios::dec | std::ios::scientific);
    observer.precision(16);
    
    
    // loop over observer time and compute Kirchhoff integral
    while (to < tf)
    {
        std::cout << "observer time to = " << to << std::endl;

        auto tauq = compute_emission_time(xq, xo, to, c);

        observer << to << "\t" << compute_kirchhoff_integral(rq, cosq, tauq, theta, t, p, pr, R, c) << std::endl;

        to += dt;
    }

    return 0;
}
