#pragma once

namespace simrail::constants {
    // Numerical tolerances
    constexpr double FLOATING_POINT_TOL = 1e-15;
    constexpr double INTEGRATION_EPSILON = 1e-20;
    constexpr double BISECTION_TOL_STRICT = 1e-10;
    constexpr double CONTACT_TOL = 1e-7;

    // Taylor expansion threshold for beta - sin(beta)
    constexpr double TAYLOR_APPROX_THRESHOLD = 1e-4;

    // Trigonometric safeguards
    constexpr double SINE_HALF_ALPHA_MIN = 1e-15;

    // Penetration correction factor
    constexpr double PENETRATION_CORRECTION = 0.55;

    // Gamma coefficient flooring (prevents force explosion)
    constexpr double GAMMA_MIN = 0.01;

    // Constants for the coefficient function Γ(Θ)
    // resembling the implemented one in SIMPACK software
    constexpr double EQU_EL_POWER_CORRECTION = 1.026600611713163413;
    constexpr double EQU_EL_DIV_CORRECTION = 1.397839912270349316;

    // Bisection iterations
    constexpr int MAX_BISECTION_ITERATIONS = 100;

    // Discretization count for contact patches
    constexpr int DEFAULT_DISCRETIZATION = 58;

    // Failure sentinel (force value that steers bisection upward)
    constexpr double SOLVER_FAILURE_FORCE = 1e10;

    // Bisection tolerance variants
    constexpr double BISECTION_TOL_DEFAULT = 1e-12;

    // Interpolation validation
    constexpr double MIN_SLOPE_DELTA_CHECK = 1e-10;
}
