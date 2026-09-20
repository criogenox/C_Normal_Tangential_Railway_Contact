#pragma once

#include <vector>

#include "Types.h"
#include "Constants.h"

namespace simrail::contact {
    /**
     * Equivalent Elastic Contact Solver
     *
     * Implements the equivalent elastic contact method for wheel-rail interaction.
     * Thread-safe for concurrent solveContact() calls (stateless after construction).
     */
    class EquivalentContactSolver {
    public:
        /**
         * @param discretization Number of points for contact patch re-discretization
         */
        explicit EquivalentContactSolver(const int discretization = constants::DEFAULT_DISCRETIZATION)
            : m_discretization(discretization) {
        }

        /**
        * Solve contact problem for given geometry and material properties.
        *
        * @param local_x       Local x-coordinates of contact patch [m]
        * @param wheel_z       Wheel profile Z-values [m]
        * @param rail_z        Rail profile Z-values [m]
        * @param radius        Local rolling radius [m]
        * @param material      Material properties
        * @param rotation_angle Contact patch rotation angle [rad]
        * @return              Contact result (may indicate no contact found)
        */
        static ContactResult solveContact(
            const std::vector<double> &local_x,
            const std::vector<double> &wheel_z,
            const std::vector<double> &rail_z,
            double radius,
            const MaterialProperties &material,
            double rotation_angle
        );

    private:
        int m_discretization;

        /// Find residual for beta root-finding (Eq. 11)
        static double findBetaResidual(double beta, double area, double width);

        /// Gamma coefficient for Hertzian correction (Eq. 21-22)
        static double gammaCoefficient(double theta_deg);
    };
} // namespace simrail::contact
