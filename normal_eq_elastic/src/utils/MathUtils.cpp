#include "MathUtils.h"

#include <algorithm>

namespace simrail::math {
    double integrateTrapezoid(const std::vector<double> &x, const std::vector<double> &y) {
        // Validation: Ensure non-degenerate state context, matching lengths
        if (x.size() < 2 || x.size() != y.size()) [[unlikely]] {
            return 0.0;
        }

        double sum = 0.0;
        for (size_t i = 1; i < x.size(); ++i) {
            const double dx = x[i] - x[i - 1];
            const double avg_y = (y[i] + y[i - 1]) * 0.5;

            // Skip negligible intervals (numerical stability)
            if (std::abs(dx) < constants::FLOATING_POINT_TOL) {
                continue;
            }

            sum += dx * avg_y;
        }
        return sum;
    }

    std::optional<double> findRootBisect(
        const std::function<double(double)> &func,
        double a, double b,
        const double tol,
        const int max_iter) {
        // Validate bounds
        if (a >= b) {
            return std::nullopt;
        }

        const double fa = func(a);
        const double fb = func(b);

        // No sign change = no guaranteed root in interval
        if (fa * fb > 0.0) return std::nullopt;

        // Early exit if one endpoint is already at root
        if (std::abs(fa) < tol) {
            return a;
        }
        if (std::abs(fb) < tol) {
            return b;
        }

        for (int i = 0; i < max_iter; ++i) {
            double mid = (a + b) / 2.0;
            const double fm = func(mid);

            // Convergence check: function value AND interval width
            if (std::abs(fm) < tol || (b - a) / 2.0 < tol)
                return mid;

            // Bisection step
            if (fa * fm < 0.0) {
                b = mid;
            } else {
                a = mid;
            }
        }

        // Returned best estimate after max iterations
        return (a + b) / 2.0;
    }
} // namespace simrail::math
