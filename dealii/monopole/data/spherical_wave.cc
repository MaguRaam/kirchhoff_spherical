#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

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

int main()
{
  // spherical wave object:
  double c0 = 250.; // wave speed
  double f0 = 100.; // frequency
  SphericalWave wave(c0, f0);

  // semicircle arc radius:
  double R = 1.0;

  // no of cells along theta direction:
  const int Ntheta = 100;

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

    //loop over theta:
    for (int i = 0; i < Ntheta; ++i)
      pfile << wave.p(R, t) << "\t" << wave.dpdr(R, t) << "\n";


    //update time and step:
    t += dt;
    step++;
  }

  return 0;
}
