#pragma once
#include <string>
#include <memory>

/**
 * @struct KalkerResultCubic
 * @brief Output structure containing the normalized contact forces and spin moment.
 *
 * Deliberately distinct from the linear library's KalkerResult so both headers can
 * be included side by side (e.g. for A/B comparison in a test driver).
 */
struct KalkerResultCubic {
    double Fx; ///< Normalized longitudinal force
    double Fy; ///< Normalized lateral force
    double Mz; ///< Normalized spin moment (aM)
};

/**
 * @class KalkerTableCubic
 * @brief Same TABCON book-of-tables lookup as KalkerTable, but using
 *        monotonicity-preserving cubic Hermite (PCHIP / Fritsch-Carlson)
 *        interpolation instead of multilinear interpolation.
 *
 * The file format, grids, and symmetry handling are identical to KalkerTable.
 * The only difference is the interpolation kernel: each of the four dimensions is
 * interpolated with a local 4-node monotone cubic, applied as a successive
 * (tensor-product) reduction phi -> uy -> ux -> aob. This captures the curvature
 * of the saturation "knee" that a straight-line interpolant misses, while the
 * monotone slope limiter guarantees the interpolant never overshoots the stored
 * data (so |F*| stays within the friction bound and the response stays monotone).
 */
class KalkerTableCubic {
public:
    KalkerTableCubic();

    ~KalkerTableCubic();

    /**
     * @brief Loads the TABCON data from a file.
     * @param filename Path to tabcon.dat.
     * @return true if loading was successful, false otherwise.
     */
    [[nodiscard]] bool load(const std::string &filename) const;

    /**
     * @brief Performs a 4D lookup with monotone cubic interpolation.
     * @param aob Semi-axis ratio (a/b), signed creepages handled internally.
     * @param ux Scaled longitudinal creepage (xi).
     * @param uy Scaled lateral creepage (eta).
     * @param phi Scaled spin (psi).
     * @return KalkerResultCubic containing interpolated Fx, Fy, and Mz.
     */
    [[nodiscard]] KalkerResultCubic lookup(double aob, double ux, double uy, double phi) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};
