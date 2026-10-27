// L/C test driver for the monotone-cubic TABCON interpolator.
//
// Loads tabcon.dat into both the multilinear library (KalkerTable) and the
// monotone-cubic library (KalkerTableCubic), runs the reference cases through
// the identical normalization, and prints linear vs cubic vs expected physical
// forces with error percentages.
//
// This file is self-contained: it re-derives the Kalker coefficients (linrol) and
// the creepage normalization locally, so it does not depend on other files.

#include "kalker_table.h"        // linear library
#include "kalker_table_cubic.h"  // cubic library

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

// ---------------------------------------------------------------------------
// Kalker linear-theory creepage coefficients (port of linrol routine).
// ---------------------------------------------------------------------------
static void linrol(const double nu, const double aob, double &cx, double &cy, double &cz) {
    const double pi = 4.0 * std::atan(1.0);

    static constexpr double cij[2][10][3][3] = {
        {
            // icase 0: a < b
            {{2.51, 3.31, 4.85}, {2.51, 2.52, 2.53}, {0.334, 0.473, 0.731}},
            {{2.59, 3.37, 4.81}, {2.59, 2.63, 2.66}, {0.483, 0.603, 0.809}},
            {{2.68, 3.44, 4.80}, {2.68, 2.75, 2.81}, {0.607, 0.715, 0.889}},
            {{2.78, 3.53, 4.82}, {2.78, 2.88, 2.98}, {0.720, 0.823, 0.977}},
            {{2.88, 3.62, 4.83}, {2.88, 3.01, 3.14}, {0.827, 0.929, 1.07}},
            {{2.98, 3.72, 4.91}, {2.98, 3.14, 3.31}, {0.930, 1.03, 1.18}},
            {{3.09, 3.81, 4.97}, {3.09, 3.28, 3.48}, {1.03, 1.14, 1.29}},
            {{3.19, 3.91, 5.05}, {3.19, 3.41, 3.65}, {1.13, 1.15, 1.40}},
            {{3.29, 4.01, 5.12}, {3.29, 3.54, 3.82}, {1.23, 1.36, 1.51}},
            {{3.40, 4.12, 5.20}, {3.40, 3.67, 3.98}, {1.33, 1.47, 1.63}},
        },
        {
            // icase 1: a > b
            {{10.7, 11.7, 12.9}, {10.7, 12.8, 16.0}, {12.2, 14.6, 18.0}},
            {{6.96, 7.78, 8.82}, {6.96, 8.14, 9.79}, {5.72, 6.63, 7.89}},
            {{5.57, 6.34, 7.34}, {5.57, 6.40, 7.51}, {3.79, 4.32, 5.01}},
            {{4.84, 5.57, 6.57}, {4.84, 5.48, 6.31}, {2.88, 3.24, 3.70}},
            {{4.37, 5.10, 6.11}, {4.37, 4.90, 5.56}, {2.35, 2.62, 2.96}},
            {{4.06, 4.78, 5.80}, {4.06, 4.50, 5.04}, {2.01, 2.23, 2.50}},
            {{3.82, 4.54, 5.58}, {3.82, 4.21, 4.67}, {1.76, 1.95, 2.18}},
            {{3.65, 4.36, 5.42}, {3.65, 3.99, 4.39}, {1.58, 1.75, 1.94}},
            {{3.51, 4.22, 5.30}, {3.51, 3.81, 4.16}, {1.44, 1.59, 1.77}},
            {{3.40, 4.12, 5.20}, {3.40, 3.67, 3.98}, {1.33, 1.47, 1.63}},
        },
    };

    int icase;
    double g;
    if (aob <= 1.0) {
        g = aob;
        icase = 0;
    } else {
        g = 1.0 / aob;
        icase = 1;
    }

    if (icase == 0 && g < 0.101) {
        // First asymptotes, g -> 0, a < b.
        cx = pi * pi / 4.0 / (1.0 - nu);
        cy = pi * pi / 4.0;
        cz = pi * std::sqrt(g) / 3.0 / (1.0 - nu) * (1.0 + nu * (std::log(16.0 / g) - 5.0));
    } else if (icase == 1 && g < 0.101) {
        // Last asymptotes, g -> 0, a > b.
        const double al = std::log(16.0 / (g * g));
        cx = 2.0 * pi / (al - 2.0 * nu) / g * (1.0 + (3.0 - std::log(4.0)) / (al - 2.0 * nu));
        cy = 2.0 * pi / g * (1.0 + (1.0 - nu) * (3.0 - std::log(4.0)) /
                             ((1.0 - nu) * al + 2.0 * nu)) /
             ((1.0 - nu) * al + 2.0 * nu);
        cz = 2.0 * pi / 3.0 / std::pow(g, 1.5) / ((1.0 - nu) * al - 2.0 + 4.0 * nu);
    } else {
        // In range of the table, 0.101 < g <= 1.0.
        // irow == floor(10*g - eps); interpolation weight fac == rem(10*g).
        int irow = static_cast<int>(10.0 * g - 0.005) - 1;
        if (irow < 0) irow = 0;
        if (irow > 8) irow = 8;
        const double fac = 10.0 * g - static_cast<double>(irow + 1);

        // Interpolate in g (row) for each nu and each coefficient.
        // NOTE: cij is laid out [icoef][inu], so index it that way.
        double c[3][3]; // c[inu][icoef]
        for (int inu = 0; inu < 3; ++inu)
            for (int icoef = 0; icoef < 3; ++icoef)
                c[inu][icoef] = fac * cij[icase][irow + 1][icoef][inu] +
                                (1.0 - fac) * cij[icase][irow][icoef][inu];

        // Interpolate to the actual nu (parabola through nu = 0, 0.25, 0.5 on 1/C).
        double cc[3];
        for (int icoef = 0; icoef < 3; ++icoef)
            cc[icoef] = (nu - 0.25) * (nu - 0.50) * 8.0 / c[0][icoef]
                        - (nu - 0.50) * nu * 16.0 / c[1][icoef]
                        + nu * (nu - 0.25) * 8.0 / c[2][icoef];

        cx = 1.0 / cc[0];
        cy = 1.0 / cc[1];
        cz = 1.0 / cc[2];
    }
}

namespace {
    // ---------------------------------------------------------------------------
    struct Case {
        // Input Physical Parameters
        double G; ///< Shear modulus (Pa)
        double nu; ///< Poisson's ratio
        double mu; ///< Friction coefficient
        double N; ///< Normal force (N)
        double a, b; ///< Contact ellipse semi-axes (m)

        // Creepages
        double gamma_x; ///< Physical longitudinal creepage (-)
        double gamma_y; ///< Physical lateral creepage (-)
        double omega_z; ///< Physical spin creepage (rad/mm)
    };

    struct Ref {
        double Fx, Fy, Mz;
    };

    struct Physical {
        double Fx, Fy, Mz;
    };
}

static double errpct(const double got, const double exp) {
    return (std::abs(exp) < 1e-12) ? 0.0 : 100.0 * (got - exp) / exp;
}

int main() {
    const KalkerTable lin;
    const KalkerTableCubic cub;

    if (!lin.load("tabcon.dat") || !cub.load("tabcon.dat")) {
        std::cerr << "Failed to load tabcon.dat. Ensure it is in the working directory." << std::endl;
        return 1;
    }
    std::cout << "Loaded tabcon.dat into both linear and cubic interpolators.\n" << std::endl;

    const Case cases[] = {
        {82e9, 0.28, 0.3, 90e3, 0.007257, 0.005155, 0.002, 0.003, 0.0001},
        {82e9, 0.28, 0.3, 90e3, 0.007257, 0.005155, 0.03, 0.02, 0.001},
        {82e9, 0.28, 0.3, 45e3, 0.007257, 0.005155, 0.003, 0.002, 0.0001},
        {82e9, 0.28, 0.3, 90e3, 0.007257, 0.005155, 0.003, 0.002, 0.0001},
    };

    const Ref refs[] = {
        {-14350.0, -22320.0, -2.51},
        {-22270.0, -15110.0, -6.393},
        {-10960.0, -7827.0, -2.605},
        {-22058.0, -15190.0, -4.022},
    };

    std::cout << std::fixed;
    for (int i = 0; i < 4; ++i) {
        const Case &c = cases[i];

        // --- PRE-PROCESSOR: Normalization ---
        // The TABCON xi, eta, psi are the Kalker-linear-theory scaled creepages
        //
        //   xi  = gamma_x * G*c^2*cx / (3*mu*N)
        //   eta = gamma_y * G*c^2*cy / (3*mu*N)
        //   psi = omega_z * G*c^3*cz / (mu*N)
        //
        // where cx, cy, cz = linrol(nu, a/b) are Kalker's Cij coefficients.
        const double L = std::sqrt(c.a * c.b); // c: characteristic length
        const double c2 = c.a * c.b; // c^2
        const double c3 = c2 * L; // c^3 = (a*b)^1.5
        const double muN = c.mu * c.N;
        // omega_z in test data is given in rad/mm -> convert to rad/m.
        const double omega_z_m = c.omega_z * 1000.0; // rad/m

        // Kalker coefficients at the PHYSICAL ratio (linrol handles a/b > 1 itself).
        double cx, cy, cz;
        linrol(c.nu, c.a / c.b, cx, cy, cz);

        // Scaled creepages = table axes. The leading minus mirrors the Fortran version
        // (cksi = -xi*..., etc.): a positive creepage produces a force that opposes
        // it, so the table is indexed with the opposite sign. lookup() then carries
        // the sign back onto Fx/Fy/Mz.
        const double xi = -c.gamma_x * c.G * c2 * cx / (3.0 * muN);
        const double eta = -c.gamma_y * c.G * c2 * cy / (3.0 * muN);
        const double psi = -omega_z_m * c.G * c3 * cz / muN;
        // The reader stores both a<b and a>b (and both phi signs), so we pass the
        // physical ratio and signed creepages directly -- no axis swap needed.
        const double aob = c.a / c.b;

        // --- Lookup with both kernels ---
        const auto rl = lin.lookup(aob, xi, eta, psi);
        const auto rc = cub.lookup(aob, xi, eta, psi);

        // --- POST-PROCESSOR: Denormalization ---
        // Physical Force = Normalized Result * mu * N
        // Physical Moment = Normalized Result * mu * N * L
        const Physical pl{.Fx = rl.Fx * muN, .Fy = rl.Fy * muN, .Mz = rl.Mz * muN * L};
        const Physical pc{.Fx = rc.Fx * muN, .Fy = rc.Fy * muN, .Mz = rc.Mz * muN * L};
        const Ref &e = refs[i];

        std::cout << "  ============================== CASE " << (i + 1) << " ==============================\n";
        std::cout << std::setprecision(4)
                << "  inputs: a/b=" << aob << " xi=" << xi << " eta=" << eta << " psi=" << psi << "\n";
        std::cout << std::setprecision(1);
        std::cout << "            " << std::setw(12) << "linear" << std::setw(12) << "cubic"
                << std::setw(12) << "expected" << std::setw(11) << "lin err%" << std::setw(11) << "cub err%" << "\n";
        std::cout << "  Fx [N]  " << std::setw(12) << pl.Fx << std::setw(12) << pc.Fx << std::setw(12) << e.Fx
                << std::setw(10) << errpct(pl.Fx, e.Fx) << "%" << std::setw(10) << errpct(pc.Fx, e.Fx) << "%\n";
        std::cout << "  Fy [N]  " << std::setw(12) << pl.Fy << std::setw(12) << pc.Fy << std::setw(12) << e.Fy
                << std::setw(10) << errpct(pl.Fy, e.Fy) << "%" << std::setw(10) << errpct(pc.Fy, e.Fy) << "%\n";
        std::cout << std::setprecision(3);
        std::cout << "  Mz [Nm] " << std::setw(12) << pl.Mz << std::setw(12) << pc.Mz << std::setw(12) << e.Mz
                << std::setprecision(1)
                << std::setw(10) << errpct(pl.Mz, e.Mz) << "%" << std::setw(10) << errpct(pc.Mz, e.Mz) << "%\n\n";
    }

    return 0;
}
