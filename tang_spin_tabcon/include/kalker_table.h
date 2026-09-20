#pragma once
#include <string>
#include <memory>

/**
 * @struct KalkerResult
 * @brief Output structure containing the normalized contact forces and spin moment.
 */
struct KalkerResult {
    double Fx; ///< Normalized longitudinal force
    double Fy; ///< Normalized lateral force
    double Mz; ///< Normalized spin moment (aM)
};

/**
 * @class KalkerTable
 * @brief Replicates the USETAB routine for multidimensional linear interpolation
 *        of Kalker's Book of Tables (TABCON).
 *
 * This class handles the loading of the TABCON data file and performs
 * quadrilinear interpolation across four dimensions:
 * 1. a/b ratio (semi-axes of the contact ellipse)
 * 2. phi (normalized spin parameter)
 * 3. ux (normalized longitudinal creepage)
 * 4. uy (normalized lateral creepage)
 */
class KalkerTable {
public:
    KalkerTable();

    ~KalkerTable();

    /**
     * @brief Loads the TABCON data from a file.
     * @param filename path to tabcon.dat.
     * @return true if loading was successful, false otherwise.
     */
    [[nodiscard]] bool load(const std::string &filename) const;

    /**
     * @brief Performs a 4D lookup and interpolation.
     *
     * This method handles internal symmetry relations and clips inputs to the grid range.
     *
     * @param aob Semi-axis ratio (a/b).
     * @param ux Normalized longitudinal creepage.
     * @param uy Normalized lateral creepage.
     * @param phi Normalized spin parameter.
     * @return KalkerResult containing interpolated Fx, Fy, and Mz.
     */
    [[nodiscard]] KalkerResult lookup(double aob, double ux, double uy, double phi) const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl;
};
