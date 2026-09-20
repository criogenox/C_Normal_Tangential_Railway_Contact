#pragma once

#include <array>

/**
 * @struct Inputs
 * @brief Input parameters and force accumulators for the FASTSIM algorithm.
 */
struct Inputs {
    double ux; ///< Normalized longitudinal creepage
    double uy; ///< Normalized lateral creepage
    double fx; ///< Normalized spin parameter in longitudinal direction
    double fy; ///< Normalized spin parameter in lateral direction
    double mx; ///< Number of strip divisions in x-direction
    double my; ///< Number of strip divisions in y-direction
    double Tx; ///< Accumulated longitudinal tangential force
    double Ty; ///< Accumulated lateral tangential force
};

/**
 * @namespace Helper
 * @brief Utility functions for tangential traction and force calculations.
 */
namespace Helper {
    /**
     * @brief Computes the updated tangential traction along a strip.
     * @param Pold Previous traction value.
     * @param S Creepage / spin term.
     * @param x_prev Previous x-coordinate.
     * @param x_curr Current x-coordinate.
     * @return Updated tangential traction.
     */
    inline double pressure(const double Pold, const double S, const double x_prev, const double x_curr) {
        return Pold - S * (x_prev - x_curr);
    }

    /**
     * @brief Accumulates tangential force contribution from a contact element.
     * @param Told Previous accumulated force.
     * @param AR Area of the contact element (dx * dy).
     * @param P Tangential traction at the element.
     * @return Updated accumulated force.
     */
    inline double force(const double Told, const double AR, const double P) {
        return Told + AR * P;
    }
}

/**
 * @class Subroutine
 * @brief Base class for FASTSIM stripwise tangential traction integration.
 */
class Subroutine {
protected:
    Inputs args_; ///< Algorithm input parameters and force accumulators

public:
    explicit Subroutine(const Inputs &args);

    /**
     * @brief Integrates tangential traction across a single strip at lateral position y.
     * @param dy Width of the strip in y-direction.
     * @param y Lateral coordinate of the strip center.
     */
    void traction(double dy, double y);
};

/**
 * @class Fastsim
 * @brief Implements Kalker's FASTSIM algorithm for the simplified theory of rolling contact.
 *
 * This class computes the normalized tangential forces (longitudinal and lateral)
 * acting on an elliptical contact patch by integrating tractions across contact strips.
 */
class Fastsim : public Subroutine {
public:
    /**
     * @throws std::invalid_argument if mx or my < 1 or if mx or my are not integers.
     */
    explicit Fastsim(const Inputs &args);

    ~Fastsim() = default;

    Fastsim(const Fastsim &other) = default;

    Fastsim(Fastsim &&other) noexcept = default;

    std::array<double, 2> T{}; ///< Output normalized tangential forces {Tx, Ty}

    /**
     * @brief Executes the FASTSIM integration over the contact ellipse.
     *
     * Iterates over contact ellipse strips using adaptive or uniform discretization
     * depending on the creepage regime and tolerance.
     *
     * @param TOL Tolerance threshold for adaptive strip refinement.
     */
    void creep(double TOL);
};
