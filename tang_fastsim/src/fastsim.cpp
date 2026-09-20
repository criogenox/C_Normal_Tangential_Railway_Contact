#include <cmath>
#include <numbers>
#include <stdexcept>

#include "fastsim.h"

using namespace Helper;

Subroutine::Subroutine(const Inputs &args) : args_(args) {
}

/* ==========================================================================
   Subroutine::traction
   --------------------------------------------------------------------------
   Integrates tangential traction over one y-strip of the elliptical contact
   patch (Kalker's simplified theory, 1982).
   - Traverses the strip from leading edge (x = +x0) to trailing edge (x = -x0).
   - Accumulates force contributions in args_.Tx / args_.Ty.
   ========================================================================== */
void Subroutine::traction(const double dy, const double y) {
    double Px = 0.0, Py = 0.0;

    const double y2 = y * y;
    double x_prev = std::sqrt(1.0 - y2); // semi-axis length at y
    const double dx = 2.0 * x_prev / args_.mx;
    const double area_element = dx * dy; // ΔA = dx * dy
    const double B = x_prev; // trailing edge = -B

    double x = x_prev - dx / 2.0; // start at leading-edge center
    const double Sx = args_.ux - y * args_.fx; // longitudinal creep + spin term

    while (x > -B) // traverse strip
    {
        const double Sy = args_.uy + args_.fy * (x_prev + x) / 2.0;

        Px = pressure(Px, Sx, x_prev, x);
        Py = pressure(Py, Sy, x_prev, x);

        // Check friction limit (normalized Hertz pressure = sqrt(1-y²-x²))
        const double denom = 1.0 - y2 - x * x;
        const double norm = std::hypot(Px, Py);
        /* Improves numerical stability ↑↑
        compared to sqrt(Px*Px + Py*Py) */
        if (norm / denom > 1.0) {
            const double scale = denom / norm;
            Px *= scale;
            Py *= scale;
        }

        // Accumulate force contribution
        args_.Tx = force(args_.Tx, area_element, Px);
        args_.Ty = force(args_.Ty, area_element, Py);

        // Data for next segment
        x_prev = x;
        x -= dx;
    }
}

Fastsim::Fastsim(const Inputs &args) : Subroutine(args) {
    if (args_.mx < 1.0 || args_.my < 1.0) {
        throw std::invalid_argument("mx and my must be >= 1");
    }
    if (std::floor(args_.mx) != args_.mx ||
        std::floor(args_.my) != args_.my) {
        throw std::invalid_argument("mx and my must be integers");
    }
}

/* ==========================================================================
   Fastsim::creep
   --------------------------------------------------------------------------
   Main integration driver:
   - Loops over upper (s = +1) and lower (s = -1) half of the contact ellipse.
   - Two regimes:
       1. Small-creep regime (adaptive strip width dy)
       2. Large-creep regime (uniform strip width, single pass)

   Added guards:
   - Near-singular ratio (|ux| ≈ |fx| or fx ≈ 0) → fall back to uniform strips.
     Prevents dy→0 infinite loops and NaN propagation.
   - Refinement cascade capped at MAX_REFINEMENTS per half to protect
     against pathological TOL values.
   ========================================================================== */
void Fastsim::creep(const double TOL) {
    constexpr double RATIO_TOL = 1e-6;
    const bool large_creep =
            std::abs(args_.ux) > std::abs(args_.fx) ||
            std::abs(args_.uy) > std::abs(args_.fy);
    const bool near_singular =
            std::abs(args_.fx) < RATIO_TOL ||
            std::abs(std::abs(args_.ux) - std::abs(args_.fx)) < RATIO_TOL;

    double dy = 0.0;
    double y_sign = 1.0; // s in original paper: +1 / -1

    while (y_sign > -2.0) // process upper then lower half
    {
        bool use_adaptive_strips;
        double y_min_for_this_half;

        if (large_creep || near_singular) {
            // Uniform strips, single pass over full height
            use_adaptive_strips = false;
            dy = 2.0 / args_.my;
            y_sign = -1.0; // force lower half and exit
            y_min_for_this_half = -1.0;
        } else {
            // Adaptive dy based on creep boundary
            use_adaptive_strips = true;
            const double ratio = args_.ux / args_.fx * y_sign;
            dy = (1.0 - ratio) /
                 (std::floor((1.0 - ratio) * args_.my / 2.0) + 1.0);
            y_min_for_this_half = ratio + dy;
        }

        constexpr int MAX_REFINEMENTS = 20;
        int refinements = 0;

        bool continue_integration = true;
        double y_strip_center = (1.0 + dy / 2.0) * y_sign;

        while (continue_integration) {
            y_strip_center -= dy * y_sign;

            if (y_strip_center * y_sign < y_min_for_this_half && !use_adaptive_strips) {
                continue_integration = false;
            } else if (y_strip_center * y_sign < y_min_for_this_half) {
                // Adaptive refinement when tolerance is tight
                if (TOL < (dy / 2.0) && refinements < MAX_REFINEMENTS) {
                    ++refinements;
                    dy /= 2.0;
                    y_strip_center += (dy / 2.0) * y_sign;
                    traction(dy, y_strip_center);
                    continue_integration = true; // repeat with finer strip
                } else {
                    traction(dy, y_strip_center);
                    continue_integration = false;
                }
            } else {
                traction(dy, y_strip_center);
            }
        }

        y_sign -= 2.0; // switch to lower half
    }

    // Final scaling to physical tangential forces (Kalker's normalization)
    T = {
        2.0 * args_.Tx / std::numbers::pi,
        2.0 * args_.Ty / std::numbers::pi
    };
}
