// Read pressure data on a spherical surface and compute the far-field pressure:

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

// create grid:
inline std::vector<double> create_grid(double xmin, double xmax, int ncells)
{
    std::vector<double> x(ncells);
    double dx = (xmax - xmin) / static_cast<double>(ncells);

    for (int i = 0; i < ncells; ++i)
        x[i] = xmin + (i + 0.5) * dx;

    return x;
}

// compute quadrature points on sphere:
auto create_sphere(const double R, const std::vector<double> &theta, const std::vector<double> &phi)
{

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

// write points on sphere in vtk format
void write_sphere(const boost::multi_array<Point, 2> &xq)
{
    std::ofstream sphere("sphere.vtk");

    sphere << "# vtk DataFile Version 2.0\n";
    sphere << "vtk output\n";
    sphere << "ASCII\n";
    sphere << "DTATASET POLYDATA\n";
    sphere << "POINTS " << xq.num_elements() << " float\n";

    for (auto iter = xq.data(); iter != xq.data() + xq.num_elements(); ++iter)
        sphere << iter->x << " " << iter->y << " " << iter->z << "\n";
}

// load pressure data:
auto load_pressure(const std::string filename)
{
    std::vector<double> t, p, pr;
    double col1, col2, col3;
    std::ifstream file(filename);

    while (file >> col1 >> col2 >> col3)
    {
        t.push_back(col1);
        p.push_back(col2);
        pr.push_back(col3);
    }

    return std::tuple(std::move(t), std::move(p), std::move(pr));
}

// compute distance and angle between quadrature point/normal and observer point:
auto compute_distance_angle(const boost::multi_array<Point, 2> &xq, Point xo)
{

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

// 3 point Lagrange interpolation at emission time:
inline auto interpolate_pressure(double tau, const std::vector<double> &t, const std::vector<double>& p, const std::vector<double>& pr)
{

    //we assume pressure at negative emission time is 0: 
    if (tau <= 0.)
        return std::tuple(0.0, 0.0, 0.0);

    // get the time stencil for interpolation:
    int i = 0;

    while (tau > t[i])
        i++;

    
    // interpolate pressure and its derivative at emission time using 3pt Lagrange polynomial:
    double p_tau = lagrange_interpolation<3>({t[i - 1], t[i], t[i + 1]}, {p[i-1], p[i], p[i + 1]}, tau);
    double pr_tau = lagrange_interpolation<3>({t[i - 1], t[i], t[i + 1]}, {pr[i-1], pr[i], pr[i + 1]}, tau);
    double pt_tau = lagrange_derivative<3>({t[i - 1], t[i], t[i + 1]}, {p[i-1], p[i], p[i + 1]} , tau);

    return std::tuple(std::move(p_tau), std::move(pr_tau), std::move(pt_tau));
}

// compute kirchhoff integral:
double compute_kirchhoff_integral(const boost::multi_array<double, 2> &rq, const boost::multi_array<double, 2> &cosq, const boost::multi_array<double, 2> &tauq, const std::vector<double> &t, const std::vector<double> &p, const std::vector<double> &pr, double R, double c, const std::vector<double> &theta)
{

    double integral = 0.0;

    // no of cells:
    int n_theta = rq.shape()[0];
    int n_phi = rq.shape()[1];

    // cell size:
    const double dtheta = M_PI / static_cast<double>(n_theta);
    const double dphi = 2. * M_PI / static_cast<double>(n_phi);

    double cinv = 1. / c;

    // loop over theta and phi grid points:
    for (int i = 0; i < n_theta; ++i)
    {
        for (int j = 0; j < n_phi; ++j)
        {

            auto [p_tau, pr_tau, pt_tau] = interpolate_pressure(tauq[i][j], t, p, pr);

            integral += (0.25 / M_PI) * ((cinv * pt_tau * cosq[i][j] - pr_tau) / rq[i][j] + (p_tau * cosq[i][j]) / (rq[i][j] * rq[i][j])) * R * R * sin(theta[i]) * dtheta * dphi;
        }
    }

    return integral;
}

int main()
{
    // wave speed in the water medium:
    const double c = 163.82917933017916;

    // radius of sphere:
    const double R = 5.001;

    // observer point:
    const Point xo{10.001, 0.0, 0.0};

    // observer time:
    double to = 0.0;
    const double tf = 4.0;
    const double dt = 0.1;

    // no of cells along theta and phi direction:
    const int n_theta = 5000;
    const int n_phi = 2 * n_theta;

    // data filename:
    std::string filename = "data/R=5.001000.dat";

    // compute cell size:
    double h = M_PI / static_cast<double>(n_theta);

    // create theta and phi grid points:
    auto theta = create_grid(0, M_PI, n_theta);
    auto phi = create_grid(0, 2. * M_PI, n_phi);

    // compute quadrature points on sphere:
    std::cout << "create sphere \n";
    auto xq = create_sphere(R, theta, phi);

    // load pressure and time data:
    std::cout << "load pressure data \n";
    auto [t, p, pr] = load_pressure(filename);

    // compute distance and angle between quadrature point/normal and observer point:
    std::cout << "compute distance and angle \n";
    auto [rq, cosq] = compute_distance_angle(xq, xo);
   
    // write Kirchhoff data:
    std::ofstream observer("p_kirchhoff.dat", std::ios::out);
    observer.flags(std::ios::dec | std::ios::scientific);
    observer.precision(16);


    // loop over observer time and compute Kirchhoff integral

    while (to < tf){

        std::cout << "observer time to = " << to << std::endl;

        auto tauq = compute_emission_time(xq, xo, to, c);

        //compute max and min emisiion time:
        //auto [min, max] = std::minmax_element(tauq.origin(), tauq.origin() + tauq.num_elements());
        //std::cout << "emission time range [ " << *min << "," << *max << "]\n";

        observer << to << "\t" << compute_kirchhoff_integral(rq, cosq, tauq, t, p, pr, R, c, theta) << std::endl;

        to += dt;
    }

    return 0;

}




