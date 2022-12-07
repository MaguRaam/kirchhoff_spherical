// Compute Kirchhoff integral using 2 point Gaussian Quadrature and show 4th order convergence 
// Using 5pt Lagrange polynomial in time:

// boost headers:
#include <boost/multi_array.hpp>

// c++ headers:
#include <iostream>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <sstream>
#include <cassert>
#include <tuple>

// Kirchhoff Data:
struct KirchhoffData
{
    double theta, p, pr;
};

// R3 vector:
template <typename T>
struct Vector
{
    T x, y, z;
};

// add two Vectors:
template <typename T>
inline Vector<T> operator+(const Vector<T> &p1, const Vector<T> &p2)
{
    return {p1.x + p2.x, p1.y + p2.y, p1.z + p2.z};
}

// subtract two Vectors:
template <typename T>
inline Vector<T> operator-(const Vector<T> &p1, const Vector<T> &p2)
{
    return {p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
}

// dot product:
template <typename T>
inline T dot(const Vector<T> &p1, const Vector<T> &p2)
{
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

// euclidean norm of a vector:
template <typename T>
inline T norm(const Vector<T> &p)
{
    return sqrt(dot(p, p));
}

// normalize vector:
template <typename T>
inline Vector<T> normalize(const Vector<T> &p)
{
    return {p.x / norm(p), p.y / norm(p), p.z / norm(p)};
}

// cos of angle between two vectors:
template <typename T>
inline T cos(const Vector<T> &p1, const Vector<T> &p2)
{
    return dot(normalize(p1), normalize(p2));
}

// scalar multiplied by vector
template <typename T>
inline Vector<T> operator*(T alpha, const Vector<T> &p)
{
    return {p.x * alpha, p.y * alpha, p.z * alpha};
}

// print Vector:
template <typename T, typename Stream>
inline std::ostream &operator<<(Stream &os, const Vector<T> &p)
{
    os << "(" << p.x << "," << p.y << "," << p.z << ")\n";
    return os;
}

// point in 3d space:
using Point = Vector<double>;

// emission time  τ = t − |x−y|/c:
inline double emission_time(const Point &x, const Point &y, double t, double c)
{
    return t - norm(x - y) / c;
}

// spherical to cartesian:
inline Point spherical_to_cartesian(double r, double theta, double phi)
{
    return {r * cos(phi) * sin(theta), r * sin(phi) * sin(theta), r * cos(theta)};
}

// compute the ith lagrange basis evaluated at tau:
template <int npts>
inline double lagrange_basis(const std::array<double, npts> &t, double tau, int i)
{
    double value = 1.;

    // compute the ith basis:
    for (int j = 0; j < i; ++j)
        value *= (tau - t[j]) / (t[i] - t[j]);

    for (int j = i + 1; j < npts; ++j)
        value *= (tau - t[j]) / (t[i] - t[j]);

    return value;
}

// compute the jth lagrange basis derivative evaluated at tau:
template <int npts>
inline double lagrange_basis_derivative(const std::array<double, npts> &t, double tau, int j)
{
    double value = 0.;

    for (int i = 0; i < npts; i++)
    {
        if (i != j)
        {
            // product:
            double prod = 1.0;

            for (int m = 0; m < npts; ++m)
                if (m != i && m != j)
                    prod *= (tau - t[m]) / (t[j] - t[m]);

            value += (1. / (t[j] - t[i])) * prod;
        }
    }

    return value;
}

// evalue lagrange interpolation at tau:
template <int npts>
double lagrange_interpolation(const std::array<double, npts> &t, const std::array<double, npts> &p, double tau)
{
    std::array<double, npts> basis;

    for (int i = 0; i < npts; ++i)
        basis[i] = lagrange_basis<npts>(t, tau, i);

    return std::inner_product(basis.begin(), basis.end(), p.begin(), 0.0);
}

// evalue lagrange derivative at tau:
template <int npts>
inline double lagrange_derivative(const std::array<double, npts> &t, const std::array<double, npts> &p, double tau)
{
    std::array<double, npts> basis;

    for (int i = 0; i < npts; ++i)
        basis[i] = lagrange_basis_derivative<npts>(t, tau, i);

    return std::inner_product(basis.begin(), basis.end(), p.begin(), 0.0);
}

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

// compute quadrature points on sphere:
auto create_sphere(const double R, const std::vector<double> &theta, const std::vector<double> &phi)
{
	//no of quadrature points along theta and phi direction:
    int n_theta = theta.size();
    int n_phi = phi.size();

    // allocate memory:
    boost::multi_array<Point, 2> xq(boost::extents[n_theta][n_phi]);

    // loop over theta and phi grid points:
    for (int i = 0; i < n_theta; ++i)
        for (int j = 0; j < n_phi; ++j)
            xq[i][j] = spherical_to_cartesian(R, theta[i], phi[j]);

    return xq;
}

// compute distance and angle between quadrature point/normal and observer point:
auto compute_distance_angle(const boost::multi_array<Point, 2> &xq, Point xo)
{
	//no of quadrature points along theta and phi direction:
    int n_theta = xq.shape()[0];
    int n_phi = xq.shape()[1];

    // allocate memory:
    boost::multi_array<double, 2> rq(boost::extents[n_theta][n_phi]), cosq(boost::extents[n_theta][n_phi]);

    // loop over theta and phi grid points:
    for (int i = 0; i < n_theta; ++i)
    {
        for (int j = 0; j < n_phi; ++j)
        {
            rq[i][j] = norm(xo - xq[i][j]);
            cosq[i][j] = cos(xq[i][j], xo - xq[i][j]);
        }
    }

    return std::tuple(std::move(rq), std::move(cosq));
}

// compute emission time:
auto compute_emission_time(const boost::multi_array<Point, 2> &xq, Point xo, double to, double co)
{
	//no of quadrature points along theta and phi direction:
    int n_theta = xq.shape()[0];
    int n_phi = xq.shape()[1];

    // allocate memory:
    boost::multi_array<double, 2> tauq(boost::extents[n_theta][n_phi]);

    // loop over theta and phi grid points:
    for (int i = 0; i < n_theta; ++i)
        for (int j = 0; j < n_phi; ++j)
            tauq[i][j] = emission_time(xo, xq[i][j], to, co);

    return tauq;
}

// load pressure and time data:
auto load_pressure(int n_theta)
{

    // load time data:
    std::ifstream tfile("data/time.dat");
    std::istream_iterator<double> start(tfile), end;
    std::vector<double> t(start, end);

    // allocate memory:
    int nt = t.size(); // no of time steps:
    boost::multi_array<KirchhoffData, 2> data(boost::extents[nt][n_theta]);

    // load pressure data:
    std::string line;
    std::ifstream pfile("data/pressure.dat");

    int k = 0, i = 0;

    //row:
    std::vector<KirchhoffData> row(n_theta);

    // skip first line:
    std::getline(pfile, line);

    while (std::getline(pfile, line))
    {

        // if you find t then increment time step and initialize i index to 0 and skip the line:
        if (line.find("t") != std::string::npos)
        {
            //sort the vector
            std::sort(row.begin(), row.end(), [](const KirchhoffData &a, const KirchhoffData &b) { return a.theta < b.theta;});

            //copy row vector to 2d array:
            for (int l = 0; l < n_theta; ++l)
                data[k][l] = row[l];


            k++;
            i = 0;
            continue;
        }

        std::istringstream ss(line);

        ss >> row[i].theta >> row[i].p >> row[i].pr;

        i++;
    }

    // sort the pressure based on theta:

    return std::tuple(std::move(t), std::move(data));
}

// 5 point Lagrange interpolation at emission time:
inline auto interpolate_pressure(const std::vector<double> &t, double tau, int theta_index, const boost::multi_array<KirchhoffData, 2> &data)
{

    //we assume pressure at negative emission time is 0: 
    if (tau <= 0.)
        return std::tuple(0.0, 0.0, 0.0);

    // get the time stencil for interpolation:
    int i = 0;

    while (tau > t[i])
        i++;

    // interpolate pressure and its derivative at emission time using 5pt Lagrange polynomial:
    double p_tau = lagrange_interpolation<5>({t[i - 2], t[i - 1], t[i], t[i + 1], t[i + 2]}, {data[i - 2][theta_index].p, data[i - 1][theta_index].p, data[i][theta_index].p, data[i + 1][theta_index].p, data[i + 2][theta_index].p}, tau);
    double pr_tau = lagrange_interpolation<5>({t[i - 2], t[i - 1], t[i], t[i + 1], t[i + 2]}, {data[i - 2][theta_index].pr, data[i - 1][theta_index].pr, data[i][theta_index].pr, data[i + 1][theta_index].pr, data[i + 2][theta_index].pr}, tau);
    double pt_tau = lagrange_derivative<5>({t[i - 2], t[i - 1], t[i], t[i + 1], t[i + 2]}, {data[i - 2][theta_index].p, data[i - 1][theta_index].p, data[i][theta_index].p, data[i + 1][theta_index].p, data[i + 2][theta_index].p}, tau);

    return std::tuple(std::move(p_tau), std::move(pr_tau), std::move(pt_tau));
}

// compute kirchhoff integral:
double compute_kirchhoff_integral(const boost::multi_array<double, 2> &rq, const boost::multi_array<double, 2> &cosq, const boost::multi_array<double, 2> &tauq, const std::vector<double> &theta, const std::vector<double> &t,
                                  const boost::multi_array<KirchhoffData, 2> &data, double R, double c)
{

    double integral = 0.0;

   // no of quadpoints:
    int n_theta = rq.shape()[0];
    int n_phi = rq.shape()[1];

    // cell size:
    const double dtheta = M_PI / static_cast<double>(0.5*n_theta);
    const double dphi = 2. * M_PI / static_cast<double>(0.5*n_phi);

    double cinv = 1. / c;

    double w_theta = 0.5, w_phi = 0.5;

    // loop over theta and phi grid points:
    for (int i = 0; i < n_theta; ++i)
    {
        for (int j = 0; j < n_phi; ++j)
        {

            auto [p_tau, pr_tau, pt_tau] = interpolate_pressure(t, tauq[i][j], i, data);

            integral += (0.25 / M_PI) * ((cinv * pt_tau * cosq[i][j] - pr_tau) / rq[i][j] + (p_tau * cosq[i][j]) / (rq[i][j] * rq[i][j])) * R * R * sin(theta[i]) * w_theta* w_phi* dtheta * dphi;
        }
    }

    return integral;
}



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

int main(){

    // wave speed:
    const double c = 250.0;

		// exact solution:
    SphericalWave wave(c, 100.0);

    // radius of sphere:
    const double R = 1.0;

    // observer point:
    const Point xo{1.5, 0.0, 0.0};

	// observer time:
    double to = 0.04;
		
    // write convergence data:
    std::ofstream file("error.dat", std::ios::app);
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    // no of cells along theta and phi direction:
    const int n_theta = 1600;
    const int n_phi = 2 * n_theta;

    // compute cell size:
    double h = M_PI / static_cast<double>(n_theta);

    // create quadrature points along theta and phi direction:
    auto theta = create_quadrature_points(0, M_PI, n_theta);
    auto phi = create_quadrature_points(0, 2. * M_PI, n_phi);

    // compute quadrature points on sphere:
    auto xq = create_sphere(R, theta, phi);

    // load pressure and time data:
    auto [t, data] = load_pressure(2*n_theta);

    // compute distance and angle between quadrature point/normal and observer point:
    auto [rq, cosq] = compute_distance_angle(xq, xo);

    // compute emission time at each quadrature point:
    auto tauq = compute_emission_time(xq, xo, to, c);

    double error = fabs(compute_kirchhoff_integral(rq, cosq, tauq, theta, t, data, R, c) - wave.p(norm(xo), to));

    file << h << "\t" << error << std::endl;
}
