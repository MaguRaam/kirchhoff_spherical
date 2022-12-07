//generate wave data from monopole source at quadrature points located at circular arc (we assume 2pt Gauss quadratue):

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <array>

// Monopole solution:
class SphericalWave
{

public:
  // constructor:
  SphericalWave(double c0, double f0) : c0(c0), f0(f0), t0(4.0 / f0) {}

  // pressure p
  double p(double r, double t)
  {

    return source(t - r / c0) / (4.0 * M_PI * r);
  }

  // pressure derivative dp/dr:
  double dpdr(double r, double t)
  {

    return source_derivative(t - r / c0) / (-c0 * 4.0 * M_PI * r) + source(t - r / c0) * (-0.25 / (M_PI * r * r));
  }

  // pressure derivative dp/dt:
  double dpdt(double r, double t)
  {

    return source_derivative(t - r / c0) / (4.0 * M_PI * r);
  }

private:
  double source(double t) { return -2.0 * (t - t0) * f0 * f0 * exp(-1.0 * f0 * f0 * (t - t0) * (t - t0)); }
  double source_derivative(double t) { return -1.0 * f0 * f0 * f0 * f0 * (-2.0 * t + 2.0 * t0) * (2 * t - 2 * t0) * exp(-1.0 * f0 * f0 * (t - t0) * (t - t0)) - 2.0 * f0 * f0 * exp(-1.0 * f0 * f0 * (t - t0) * (t - t0)); }

  const double c0;
  const double f0;
  const double t0;
};

// create quadrature points along theta and phi direction:
std::vector<double> create_quadrature_points(double xmin, double xmax, int ncells)
{
    std::vector<double> xq(2 * ncells);
    double dx = (xmax - xmin) / static_cast<double>(ncells);

    // 2pt quadrature on reference cell [-0.5, 0.5]:
    std::array<double, 2> qpt{-0.28867513459481287, 0.28867513459481287};

    // map quadrature points on the computational cell:
    for (int i = 0; i < ncells; ++i)
    {
        // compute cell center:
        double xc = xmin + (i + 0.5) * dx;

        // compute coordinates of quadrature point:
        xq[2 * i] = xc + dx * qpt[0];
        xq[2 * i + 1] = xc + dx * qpt[1];
    }

    return xq;
}

int main()
{
    // spherical wave object:
    double c0 = 250.; // wave speed
    double f0 = 100.; // frequency
    SphericalWave wave(c0, f0);

    // semicircle arc radius:
    double R = 1.0;

    // no of cells along theta direction:
    const int Ntheta = 200;

    // get theta at quadrature points:
    const auto theta = create_quadrature_points(0.0, M_PI, Ntheta);

    // write time data:
    int step = 0;
    double t = 0.0, tf = 0.1, dt = 0.0001;

    //open time file:
    std::ofstream pfile("pressure.dat"), tfile("time.dat");
    pfile.flags(std::ios::dec | std::ios::scientific);
    tfile.flags(std::ios::dec | std::ios::scientific);
    pfile.precision(16);
    tfile.precision(16);

    // evolve in time:
    while (t < tf)
    {
        // write time data:
        pfile << "t = "<< t << "\n";
        tfile << t << "\n";

        //loop over cells:
        for (int i = 0; i < Ntheta; ++i){

            //first quadrature point:
            pfile << theta[2*i] << "\t" << wave.p(R, t) << "\t" << wave.dpdr(R, t) << "\n";

            //second quadrature point:
            pfile << theta[2*i + 1] << "\t"<< wave.p(R, t) << "\t" << wave.dpdr(R, t) << "\n";

        }

        //update time and step:
        t += dt;
        step++;
    }
}
