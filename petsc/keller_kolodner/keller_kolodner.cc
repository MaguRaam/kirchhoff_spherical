// Solving for bubble radius and far-field pressure using the Keller Kolodner model:

#include <iostream>
#include <fstream>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

using state_type = std::vector<double>;

// simulation parameters:

// medium
const double p_inf = 100;                                          // Pressure of liquid far away from the bubble
const double rho_inf = 1.0;                                        // Density of liquid
const double Gamma_inf = 4.4;                                      // Ratio of specific heats of the liquid
const double pi_inf = 6000.0;                                      // Stiffness constant of the liquid
const double c_inf = sqrt(Gamma_inf * (p_inf + pi_inf) / rho_inf); // Speed of sound in the liquid

// bubble:
const double p0 = 1.0;    // Initial pressure inside the bubble
const double Gamma = 1.4; // Ratio of specific heats inside the bubble
const double R0 = 1.0;    // Initial radius of the bubble
const double R0dot = 0.0; // Initial velocity of the bubble interface

const double ro = 100.001000; // observer point to compute far-field pressure

const double ti = 0.0;    // initial time
const double tf = 1.5;    // final time
const double dt = 1.0e-3; // time step

// compute the ith Lagrange basis evaluated at tau:
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

// evaluate Lagrange interpolation at tau:
template <int npts>
double lagrange_interpolation(const std::array<double, npts> &t, const std::array<double, npts> &p, double tau)
{
    std::array<double, npts> basis;

    for (int i = 0; i < npts; ++i)
        basis[i] = lagrange_basis<npts>(t, tau, i);

    return std::inner_product(basis.begin(), basis.end(), p.begin(), 0.0);
}

// interpolate R, Rdot and Rddot at emission time using 3pt Lagrange polynomial:
inline auto interpolate(double tau, const std::vector<double> &t, const std::vector<double> &R, const std::vector<double> &Rdot, const std::vector<double> &Rddot)
{
    // we assume R, Rdot, Rddot at negative emission time is 0:
    if (tau <= 0.)
        return std::tuple(0.0, 0.0, 0.0);

    // get the time stencil for interpolation:
    int i = 0;

    while (tau > t[i])
        i++;

    // interpolate R, Rdot and Rddot at emission time using 3pt Lagrange polynomial:
    double R_tau = lagrange_interpolation<3>({t[i - 1], t[i], t[i + 1]}, {R[i - 1], R[i], R[i + 1]}, tau);
    double Rdot_tau = lagrange_interpolation<3>({t[i - 1], t[i], t[i + 1]}, {Rdot[i - 1], Rdot[i], Rdot[i + 1]}, tau);
    double Rddot_tau = lagrange_interpolation<3>({t[i - 1], t[i], t[i + 1]}, {Rddot[i - 1], Rddot[i], Rddot[i + 1]}, tau);

    return std::tuple(std::move(R_tau), std::move(Rdot_tau), std::move(Rddot_tau));
}

// acoustic pressure from the oscillating bubble (Note: the inputs are given at emission time):
double pressure(double R, double Rdot, double Rddot, double r)
{
    double p_B = p0 * pow(R0 / R, 3.0 * Gamma);
    double delta = (p_B - p_inf) / rho_inf;

    double f = -R * R * Rdot + ((R * R) / c_inf) * (0.5 * Rdot * Rdot + delta);
    double fprime = -R * (0.5 * Rdot * Rdot + delta);

    // assuming c-1 we get p-p0 = -rho_inf * (-(2.0 * R * Rdot * Rdot + R * R * Rddot) / r + (pow(R, 4.0) * Rdot * Rdot) / (2.0 * pow(r, 4.0)));

    return rho_inf * ((-fprime / r) - (f * f) / (2 * pow(r, 4)) - (0.5 / c_inf) * ((fprime * fprime) / (c_inf * r * r) + (2.0 * f * fprime) / (r * r * r)));
}

int main()
{

    // right handside, x' = f(x), x = (R, dR/dt):
    auto keller_kolodner = [](const state_type &x, state_type &dxdt, const double /* t */)
    {
        // x[0] -> radius of bubble (R)
        // x[1] -> velocity of bubble interface   (dR/dt)

        // compute bubble pressure and time derivative of bubble pressure:
        double p_B = p0 * pow(R0 / x[0], 3.0 * Gamma);
        double dp_Bdt = -3.0 * Gamma * p0 * pow(R0, 3 * Gamma) * pow(x[0], (-3.0 * Gamma - 1.0)) * x[1];

        double term1 = (1.0 - x[1] / c_inf) * x[0];
        double term2 = ((3. / 2.) * (1. - x[1] / (3.0 * c_inf))) * x[1] * x[1];
        double term3 = ((1.0 + x[1] / c_inf) * (p_B - p_inf) + (x[0] / c_inf) * dp_Bdt) / rho_inf;

        dxdt[0] = x[1];
        dxdt[1] = (term3 - term2) / term1;
    };

    // initial condition:
    state_type x{R0, R0dot};

    // container to store radius, velocity, acceleration and time:
    vector<double> R, Rdot, Rddot, time;

    // observer:
    auto push_back_state_and_time = [keller_kolodner, &R, &Rdot, &Rddot, &time](const state_type &x, double t)
    {
        state_type dxdt(2);
        keller_kolodner(x, dxdt, t);

        R.push_back(x[0]);
        Rdot.push_back(x[1]);
        Rddot.push_back(dxdt[1]);
        time.push_back(t);
    };

    // ode solver type:
    using error_stepper_type = runge_kutta_cash_karp54<state_type>;

    // solve:
    size_t steps = integrate_adaptive(make_controlled(1.0e-12, 1.0e-12, error_stepper_type{}), keller_kolodner, x, ti, tf, dt, push_back_state_and_time);

    // write data:
    std::ofstream file("solution.dat");
    file.flags(std::ios::dec | std::ios::scientific);
    file.precision(16);

    // write t, R(t) and p(ro, t):
    for (size_t i = 0; i <= steps; i++)
    {
        // compute retarded time:
        double tau = time[i] - ro / c_inf;

        // interpolate R, Rdot, Rddot at emission time:
        auto [R_tau, Rdot_tau, Rddot_tau] = interpolate(tau, time, R, Rdot, Rddot);

        // compute far-field pressure:
        double p = pressure(R_tau, Rdot_tau, Rddot_tau, ro);

        file << time[i] << '\t' << R[i] << '\t' << p << '\n';
    }

    return 0;
}
