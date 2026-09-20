#pragma once

#include <cmath>
#include <functional>
#include <optional>
#include <vector>

#include "Constants.h"

namespace simrail::math {
    /**
     * Trapezoidal integration over x-y pairs.
     * Returns 0.0 if fewer than 2 points provided.
     */
    double integrateTrapezoid(const std::vector<double> &x, const std::vector<double> &y);

    /**
     * Bisection root finder — robust, no derivatives needed.
     *
     * @param func      Function to find root of (must cross zero between a and b)
     * @param a         Lower bound of search interval
     * @param b         Upper bound of search interval
     * @param tol       Convergence tolerance (absolute value of f(mid))
     * @param max_iter  Maximum iterations before aborting
     * @return          Root if found, std::nullopt if no sign change detected or converged
     */
    std::optional<double> findRootBisect(
        const std::function<double(double)> &func,
        double a, double b,
        double tol = constants::BISECTION_TOL_DEFAULT,
        int max_iter = constants::MAX_BISECTION_ITERATIONS
    );

    /**
     * Clamp value to [min_val, max_val] range.
     */
    constexpr double clamp(
        const double value,
        const double min_val,
        const double max_val) noexcept {
        return value < min_val ? min_val : (value > max_val ? max_val : value);
    }

    /**
     * Safe division with fallback for near-zero denominators.
     *
     * @param numerator       Dividend
     * @param denominator     Divisor (checked against FLOATING_POINT_TOL)
     * @param default_value   Value returned when denominator is effectively zero
     * @return                numerator / denominator, or default_value
     */
    constexpr double safeDivide(
        const double numerator,
        const double denominator,
        const double default_value = 0.0) noexcept {
        return std::abs(denominator) < constants::FLOATING_POINT_TOL ? default_value : numerator / denominator;
    }
} // namespace simrail::math
