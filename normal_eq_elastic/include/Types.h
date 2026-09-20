#pragma once

#include <vector>

namespace simrail {
    struct Point2D {
        double y;
        double z;
    };

    struct Profile {
        std::vector<double> y;
        std::vector<double> z;
        [[nodiscard]] size_t size() const noexcept { return y.size(); }
        [[nodiscard]] bool empty() const noexcept { return y.empty(); }
    };

    struct MaterialProperties {
        double elastic_modulus = 210e9;
        double shear_coefficient = 0.28;
        double friction = 0.4;

        [[nodiscard]] double effectiveModulus() const noexcept {
            return 1.0 / ((1.0 + shear_coefficient * shear_coefficient) / elastic_modulus * 2.0);
        }
    };

    struct ContactResult {
        bool contact_found = false;
        double normal_force = 0.0;
        double q_force = 0.0;
        double y_force = 0.0;
        double contact_angle = 0.0;
        double approach = 0.0;
        double semi_axis_a = 0.0;
        double semi_axis_b = 0.0;
        double centroid_y = 0.0;
        double a_b_ratio = 0.0;
        double A_hertz = 0.0;
        double B_hertz = 0.0;
        double theta_hertz = 0.0;
    };
} // namespace simrail
