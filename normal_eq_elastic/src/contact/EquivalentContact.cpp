#include "EquivalentContact.h"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "MathUtils.h"
#include "Constants.h"

namespace simrail::contact {
    ContactResult EquivalentContactSolver::solveContact(
        const std::vector<double> &local_x,
        const std::vector<double> &wheel_z,
        const std::vector<double> &rail_z,
        const double radius,
        const MaterialProperties &material,
        const double rotation_angle
    ) {
        ContactResult result;

        // Validation: matching sizes
        if (local_x.size() != wheel_z.size() || local_x.size() != rail_z.size())
            return result;

        // Validation: minimum points
        if (local_x.size() < 2) return result;

        const size_t n = local_x.size();

        // --------------------------------------------------------------------------
        // Eq. 8: Shape function S(y) = z_wheel - z_rail
        // --------------------------------------------------------------------------
        std::vector<double> shape(n);
        for (size_t i = 0; i < n; ++i)
            shape[i] = wheel_z[i] - rail_z[i];

        double width = std::abs(local_x.front() - local_x.back());

        // --------------------------------------------------------------------------
        // Eq. 9: Area of interpenetration shape
        // --------------------------------------------------------------------------
        double area = math::integrateTrapezoid(local_x, shape);
        if (area <= constants::INTEGRATION_EPSILON || width <= constants::FLOATING_POINT_TOL) return result;

        // --------------------------------------------------------------------------
        // Eq. 11: Solve for beta via bisection
        // --------------------------------------------------------------------------
        auto betaFunc = [area, width](const double beta) {
            return findBetaResidual(beta, area, width);
        };

        const auto betaOpt = math::findRootBisect(
            betaFunc,
            constants::FLOATING_POINT_TOL,
            std::numbers::pi - constants::FLOATING_POINT_TOL,
            constants::BISECTION_TOL_STRICT,
            constants::MAX_BISECTION_ITERATIONS);

        if (!betaOpt.has_value()) return result;

        const double beta = betaOpt.value();

        // --------------------------------------------------------------------------
        // Eq. 12a: Circumference radius
        // --------------------------------------------------------------------------
        const double radius_eqv = width / (2.0 * std::sin(beta / 2.0));

        // --------------------------------------------------------------------------
        // Eq. 12b: Penetration of circle segment
        // --------------------------------------------------------------------------
        const double sin_beta_4 = std::sin(beta / 4.0);
        const double approach_0 = 2.0 * radius_eqv * (sin_beta_4 * sin_beta_4); // Optimization: Avoided std::pow

        // --------------------------------------------------------------------------
        // Eq. 13: Equivalent penetration
        // --------------------------------------------------------------------------
        const double approach = approach_0 / constants::PENETRATION_CORRECTION;
        result.approach = approach;

        // --------------------------------------------------------------------------
        // Eq. 15: Leading edge xl = sqrt(2 * Rlocal * S)
        // --------------------------------------------------------------------------
        std::vector<double> xl(n);
        for (size_t i = 0; i < n; ++i)
            xl[i] = std::sqrt(2.0 * radius * std::max(shape[i], 0.0));

        // --------------------------------------------------------------------------
        // Eq. 16: Semi-axes
        // --------------------------------------------------------------------------
        result.semi_axis_a = *std::ranges::max_element(xl);
        result.semi_axis_b = width / 2.0;
        result.a_b_ratio = math::safeDivide(result.semi_axis_a, result.semi_axis_b);

        // --------------------------------------------------------------------------
        // Eq. 18: Hertz A, B, Theta
        // --------------------------------------------------------------------------
        result.A_hertz = approach / (result.semi_axis_a * result.semi_axis_a);
        result.B_hertz = approach / (result.semi_axis_b * result.semi_axis_b);

        double cosTheta = math::safeDivide(
            result.A_hertz - result.B_hertz,
            result.A_hertz + result.B_hertz, 0.0
        );
        cosTheta = math::clamp(cosTheta, -1.0, 1.0);
        result.theta_hertz = std::acos(cosTheta) * 180.0 / std::numbers::pi;

        // --------------------------------------------------------------------------
        // Eq. 21-22: Gamma(Theta) coefficient
        // --------------------------------------------------------------------------
        const double gamma = gammaCoefficient(result.theta_hertz);

        // --------------------------------------------------------------------------
        // Eq. 19: Normal force
        // --------------------------------------------------------------------------
        const double ratio = approach / gamma;
        const double ratio_cubed = ratio * ratio * ratio; // Optimization: Exponentiation via multiplication
        const double termInsideSqrt = (1.0 / (result.A_hertz + result.B_hertz)) * ratio_cubed;

        result.normal_force = (4.0 / 3.0) * material.effectiveModulus() * std::sqrt(termInsideSqrt);

        // --------------------------------------------------------------------------
        // Eq. 17: Centroid Cy = ∫(xl * y) dy / ∫xl dy
        // --------------------------------------------------------------------------
        std::vector<double> xl_times_y(n);
        for (size_t i = 0; i < n; ++i)
            xl_times_y[i] = xl[i] * local_x[i];

        const double centroidNum = math::integrateTrapezoid(local_x, xl_times_y);
        const double centroidDen = math::integrateTrapezoid(local_x, xl);
        result.centroid_y = (centroidDen > constants::INTEGRATION_EPSILON)
                                ? centroidNum / centroidDen
                                : (local_x.front() + local_x.back()) / 2.0;

        // --------------------------------------------------------------------------
        // Eq. 24: Global force components
        // --------------------------------------------------------------------------
        result.contact_angle = rotation_angle;
        result.q_force = result.normal_force * std::cos(rotation_angle);
        result.y_force = result.normal_force * std::sin(rotation_angle);
        result.contact_found = true;
        return result;
    }

    double EquivalentContactSolver::findBetaResidual(double beta, double area, double width) {
        // Eq. 10: Ac = 0.5 * (w²/4) * (β - sinβ) / sin²(β/2)
        // Simplify: Ac = (w²/8) * (β - sinβ) / sin²(β/2)
        const double sinHalfBeta = std::sin(beta / 2.0);

        if (std::abs(sinHalfBeta) < constants::SINE_HALF_ALPHA_MIN) return area; // Ac ≈ 0 when β → 0

        double beta_minus_sin_beta;
        if (beta < constants::TAYLOR_APPROX_THRESHOLD) {
            const double beta_sq = beta * beta;
            beta_minus_sin_beta = (beta_sq * beta) / 6.0;
        } else {
            beta_minus_sin_beta = beta - std::sin(beta);
        }

        const double ac = (width * width / 8.0)
                          * beta_minus_sin_beta
                          / (sinHalfBeta * sinHalfBeta);

        // Eq. 11: residual = area - Ac = 0
        return area - ac;
    }

    double EquivalentContactSolver::gammaCoefficient(double theta_deg) {
        // Eq. 21: Γ(Θ) = sqrt(1 - ((Θ-90)/90)²)^C1 / C2
        const double normalized = math::clamp((theta_deg - 90.0) / 90.0, -1.0, 1.0);
        const double term_sq = 1.0 - normalized * normalized;

        const double term = std::sqrt(term_sq);
        const double gamma = std::pow(term, constants::EQU_EL_POWER_CORRECTION) / constants::EQU_EL_DIV_CORRECTION;

        return std::max(gamma, constants::GAMMA_MIN);
    }
} // namespace simrail::contact
