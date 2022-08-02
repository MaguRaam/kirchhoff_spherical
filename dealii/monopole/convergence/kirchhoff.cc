#include "kirchhoff.h"




int main()
{

    // wave speed:
    const double c = 250.0;

    // exact solution:
    SphericalWave wave(c, 100.);

    // radius of sphere:
    const double R = 1.0;

    // observer point:
    const Point xo{3.0, 3.0, 3.0};

    // observer time:
    double to = 0.04;

    // write convergence data:
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    // no of cells along theta and phi direction:
    int n_theta = 3200;
    int n_phi = 2 * n_theta;

    // compute cell size:
    double h = M_PI / static_cast<double>(n_theta);

    // create theta and phi grid points:
    auto theta = create_grid(0, M_PI, n_theta);
    auto phi = create_grid(0, 2. * M_PI, n_phi);

    // compute quadrature points on sphere:
    auto xq = create_sphere(R, theta, phi);

    // load pressure and time data:
    auto [t, p, pr] = load_pressure(n_theta);

    // compute distance and angle between quadrature point/normal and observer point:
    auto [rq, cosq] = compute_distance_angle(xq, xo);

    // compute emission time at each quadrature point:
    auto tauq = compute_emission_time(xq, xo, to, c);

    // compute Kirchhoff integral:
    double error = fabs(compute_kirchhoff_integral(rq, cosq, tauq, theta, t, p, pr, R, c) - wave.p(norm(xo), to));

    file << h << "\t" << error << std::endl;
    

    return 0;
}
