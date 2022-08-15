#include "kirchhoff.h"

//Note there is a catch in this code emission time must be greater than 0:
//TODO add a c++ assertion to this

int main()
{

    // wave speed:
    const double c = 162.61611236282832;

    // radius of sphere:
    const double R = 0.22799999999999998;

    // observer point:
    const Point xo{0.34276, 0.0, 0.0007600000000000384};

    // observer time:
    double to = 0.00;
    const double tf = 0.007;
    const double dt = 0.0001;

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
