#pragma once

#include <vector>

#include "Types.h"

namespace simrail::geometry {
    /**
     * Piecewise Cubic Hermite Interpolating Polynomial (PCHIP).
     * Monotonicity-preserving cubic interpolation.
     *
     * Requirements:
     * - Input x values must be strictly increasing
     * - x and y vectors must have equal size >= 2
     */
    class PCHIPInterpolator {
    public:
        /**
         * @throws std::invalid_argument if x not strictly increasing or size < 2
         */
        PCHIPInterpolator(const std::vector<double> &x, const std::vector<double> &y);

        /**
         * Evaluate at single point.
         * Extrapolates linearly outside domain.
         */
        double operator()(double x) const;

        /**
         * Batch evaluation (optimized for multiple queries).
         *
         * @param x_vals Query points (unsorted, may contain duplicates)
         * @return       Interpolated y values
         */
        [[nodiscard]] std::vector<double> evaluateBatch(const std::vector<double> &x_vals) const;

    private:
        std::vector<double> m_x, m_y, m_slope, m_delta;

        void computeDerivatives();

        [[nodiscard]] double evaluateSingle(double x) const;
    };

    /**
     * Interpolate profile Z values onto target Y grid.
     *
     * @param source     Source profile (must have sorted Y)
     * @param target_y   Target Y coordinates (any order, any spacing)
     * @return           New profile with interpolated Z values
     */
    Profile interpolateProfileOntoGrid(
        const Profile &source,
        const std::vector<double> &target_y);
} // namespace simrail::geometry
