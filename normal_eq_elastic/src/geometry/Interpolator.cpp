#include "Interpolator.h"

#include <algorithm>
#include <cassert>
#include <iostream>

#include "MathUtils.h"

namespace simrail::geometry {
    PCHIPInterpolator::PCHIPInterpolator(const std::vector<double> &x, const std::vector<double> &y)
        : m_x(x), m_y(y) {
        // Validation: minimum size
        assert(m_x.size() == m_y.size() && m_x.size() >= 2);
        computeDerivatives();
    }

    void PCHIPInterpolator::computeDerivatives() {
        const size_t n = m_x.size();
        m_delta.resize(n - 1);
        m_slope.resize(n);

        // Compute segment slopes (delta_i)
        for (size_t i = 0; i < n - 1; ++i) {
            const double dx = m_x[i + 1] - m_x[i];
            const double dy = m_y[i + 1] - m_y[i];
            m_delta[i] = math::safeDivide(dy, dx, 0.0);
        }

        // Left endpoint slope
        if (std::abs(m_delta[0]) > constants::FLOATING_POINT_TOL) {
            const double denom = (m_delta[0] + (n > 2 ? m_delta[1] : m_delta[0])) * 0.5;
            if (std::abs(denom) > constants::FLOATING_POINT_TOL) {
                m_slope[0] = (m_delta[0] * m_delta[0]) / denom;

                // Clamp to avoid overshoot (PCHIP monotonicity preservation)
                if ((m_delta[0] > 0 && m_slope[0] > constants::MIN_SLOPE_DELTA_CHECK * m_delta[0]) ||
                    (m_delta[0] < 0 && m_slope[0] < constants::MIN_SLOPE_DELTA_CHECK * m_delta[0])) {
                    m_slope[0] = constants::MIN_SLOPE_DELTA_CHECK * m_delta[0];
                }
            } else {
                m_slope[0] = 0.0;
            }
        } else {
            m_slope[0] = 0.0;
        }

        // Right endpoint slope
        if (n >= 3 && std::abs(m_delta[n - 2]) > constants::FLOATING_POINT_TOL) {
            const double denom = ((m_delta[n - 3] + m_delta[n - 2]) * 0.5);
            if (std::abs(denom) > constants::FLOATING_POINT_TOL) {
                m_slope[n - 1] = (m_delta[n - 2] * m_delta[n - 2]) / denom;

                // Clamp to avoid overshoot
                if ((m_delta[n - 2] > 0 && m_slope[n - 1] > constants::MIN_SLOPE_DELTA_CHECK * m_delta[n - 2]) ||
                    (m_delta[n - 2] < 0 && m_slope[n - 1] < constants::MIN_SLOPE_DELTA_CHECK * m_delta[n - 2])) {
                    m_slope[n - 1] = constants::MIN_SLOPE_DELTA_CHECK * m_delta[n - 2];
                }
            } else {
                m_slope[n - 1] = m_delta[n - 2];
            }
        } else {
            m_slope[n - 1] = (n >= 2) ? m_delta[n - 2] : 0.0;
        }

        // Interior points
        for (size_t i = 1; i < n - 1; ++i) {
            if (m_delta[i - 1] * m_delta[i] <= 0.0) {
                // Sign change or zero => slope is zero (preserve monotonicity)
                m_slope[i] = 0.0;
            } else {
                // Harmonic mean weighted by adjacent slopes
                const double w1 = 2.0 * m_delta[i] + m_delta[i - 1];
                const double w2 = m_delta[i] + 2.0 * m_delta[i - 1];
                const double denom = (w1 / m_delta[i - 1]) + (w2 / m_delta[i]);

                if (std::abs(denom) > constants::FLOATING_POINT_TOL) {
                    m_slope[i] = (w1 + w2) / denom;
                } else {
                    m_slope[i] = 0.0;
                }
            }
        }
    }

    double PCHIPInterpolator::evaluateSingle(const double x) const {
        // Optimization: cache front/back for fast bounds checking
        const double x_min = m_x.front();
        const double x_max = m_x.back();

        // Clamp to domain boundaries (extrapolate to boundary values)
        if (x <= x_min) {
            return m_y.front();
        }
        if (x >= x_max) {
            return m_y.back();
        }

        // Binary search for interval (efficient for batch queries)
        const auto it = std::ranges::upper_bound(m_x, x);
        size_t idx = static_cast<size_t>(it - m_x.begin()) - 1;
        idx = std::clamp(idx, static_cast<size_t>(0), m_x.size() - 2);

        const double h = m_x[idx + 1] - m_x[idx];
        const double t = (x - m_x[idx]) / h;

        // Hermite basis functions
        const double h00 = (1 + 2 * t) * (1 - t) * (1 - t);
        const double h10 = t * (1 - t) * (1 - t);
        const double h01 = t * t * (3 - 2 * t);
        const double h11 = t * t * (t - 1);

        return h00 * m_y[idx]
               + h10 * m_slope[idx] * h
               + h01 * m_y[idx + 1]
               + h11 * m_slope[idx + 1] * h;
    }

    double PCHIPInterpolator::operator()(const double x) const {
        return evaluateSingle(x);
    }

    std::vector<double> PCHIPInterpolator::evaluateBatch(const std::vector<double> &x_vals) const {
        std::vector<double> result(x_vals.size());

        if (x_vals.empty()) return result;

        // Optimization: cache front/back for fast bounds checking
        const double x_min = m_x.front();
        const double x_max = m_x.back();

        for (size_t i = 0; i < x_vals.size(); ++i) {
            // Fast path for out-of-bounds (no binary search needed)
            if (const double x = x_vals[i]; x <= x_min) {
                result[i] = m_y.front();
            } else if (x >= x_max) {
                result[i] = m_y.back();
            } else {
                result[i] = evaluateSingle(x);
            }
        }
        return result;
    }

    Profile interpolateProfileOntoGrid(const Profile &source, const std::vector<double> &target_y) {
        const PCHIPInterpolator interp(source.y, source.z);

        Profile result;
        result.y = target_y;
        result.z = interp.evaluateBatch(target_y);

        return result;
    }
} // namespace simrail::geometry
