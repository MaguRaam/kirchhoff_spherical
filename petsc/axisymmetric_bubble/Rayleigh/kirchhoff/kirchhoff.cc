#include "kirchhoff.h"

int main()
{

    // wave speed:
    const double c = 51.423729930840295;

		// bubble radius:
		const double R0 = 1.0;
	
    // radius of sphere:
    const double R = 5*R0;

    // observer point:
    const Point xo{7*R0, 0.0, 0.0};

    // observer time:
    double to = 0.00;
    const double tf = 0.1;
    const double dt = 1.0;

    // no of cells along theta and phi direction:
    const int n_theta = 500;
    const int n_phi = 2 * n_theta;

    // create theta and phi grid points:
    auto theta = create_grid(0, M_PI, n_theta);
    auto phi = create_grid(0, 2. * M_PI, n_phi);

    // compute quadrature points on sphere:
    std::cout << "create sphere \n";
    auto xq = create_sphere(R, theta, phi);

    // write sphere in vtk:
    //write_sphere(xq);

    // load pressure and time data:
    std::cout << "load pressure data \n";
    auto [t, data] = load_pressure(n_theta);

    // compute distance and angle between quadrature point/normal and observer point:
    std::cout << "compute distance and angle \n";
    auto [rq, cosq] = compute_distance_angle(xq, xo);

    // write Kirchhoff data:
    std::ofstream observer("p_kirchhoff.dat", std::ios::out);
    observer.flags(std::ios::dec | std::ios::scientific);
    observer.precision(16);

    // loop over observer time and compute Kirchhoff integral
    while (to < tf)
    {
        std::cout << "compute integral at observer time to = " << to << std::endl;

        auto tauq = compute_emission_time(xq, xo, to, c);

        observer << to << "\t" << compute_kirchhoff_integral(rq, cosq, tauq, theta, t, data, R, c) << std::endl;

        to += dt;
    }

    return 0;
}
