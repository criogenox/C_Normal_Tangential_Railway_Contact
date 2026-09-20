#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <numeric>
#include <optional>
#include <utility>

#include "EquivalentContact.h"
#include "Interpolator.h"
#include "MathUtils.h"
#include "ProfileLoader.h"
#include "Types.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace simrail;

namespace {
    // ---------------------------------------------------------
    // wheel-rail side
    // ---------------------------------------------------------
    enum class Side { LEFT, RIGHT };
}

// ---------------------------------------------------------
// Resample profile: PCHIP resample to ~20000 points + negate z
// Matches Python calculate_profile() with DeltaX=0.025
// ---------------------------------------------------------
static Profile resampleProfile(const Profile &raw) {
    constexpr double DeltaX = 1.0 / 40.0;
    const size_t n_raw = raw.y.size();
    std::vector<double> idx_new;

    std::vector<double> idx_old(n_raw);
    std::iota(idx_old.begin(), idx_old.end(), 0.0);

    // 1. Calculate the maximum value x can reach based on your original condition: x < n_raw - DeltaX
    const double max_x_limit = static_cast<double>(n_raw) - DeltaX;

    // 2. Determine total steps. Added a tiny epsilon (1e-9) to protect against floating-point precision drops.
    size_t total_steps = 0;
    if (max_x_limit > 0.0) {
        total_steps = static_cast<size_t>(std::ceil(max_x_limit / DeltaX - 1e-9));
    }

    // 3. Optimize memory allocation
    idx_new.reserve(total_steps);

    // 4. Safe, clean, warning-free loop
    for (size_t i = 0; i < total_steps; ++i) {
        idx_new.push_back(static_cast<double>(i) * DeltaX);
    }

    const geometry::PCHIPInterpolator y_interp(idx_old, raw.y);
    const std::vector<double> y_new = y_interp.evaluateBatch(idx_new);

    const geometry::PCHIPInterpolator z_interp(raw.y, raw.z);
    const std::vector<double> z_new = z_interp.evaluateBatch(y_new);

    Profile result;
    result.y = y_new;
    result.z.resize(z_new.size());
    for (size_t i = 0; i < z_new.size(); ++i)
        result.z[i] = -z_new[i];

    return result;
}

// ---------------------------------------------------------
// Rotate a 2D point (y, z) by angle
// ---------------------------------------------------------
static void rotatePoint(double &y, double &z, double angle) {
    const double ca = std::cos(angle);
    const double sa = std::sin(angle);
    const double ny = ca * y - sa * z;
    const double nz = sa * y + ca * z;
    y = ny;
    z = nz;
}

// ---------------------------------------------------------
// High-Performance Linear interpolation (with stateful caching hint)
// ---------------------------------------------------------
static double linearInterp(const std::vector<double> &xp,
                           const std::vector<double> &fp,
                           const double x) {
    if (x <= xp.front()) return fp.front();
    if (x >= xp.back()) return fp.back();

    // Cache the last matched index across calls since input grids are correlated
    thread_local size_t last_idx = 0;

    // Bounds check cache safety
    if (last_idx >= xp.size() - 1) last_idx = 0;

    // Check if the point falls into the identical or adjacent neighborhood first (O(1))
    if (x >= xp[last_idx] && x <= xp[last_idx + 1]) {
        // Cache hit
    } else {
        // Fallback to binary search lookup only on cache miss
        const auto it = std::ranges::upper_bound(xp, x);
        last_idx = static_cast<size_t>(it - xp.begin()) - 1;
    }

    const double t = (x - xp[last_idx]) / (xp[last_idx + 1] - xp[last_idx]);
    return fp[last_idx] + t * (fp[last_idx + 1] - fp[last_idx]);
}

namespace {
    // ---------------------------------------------------------
    // Wheel setup — matches Python wheel class
    // ---------------------------------------------------------
    struct WheelSetup {
        std::vector<double> y_pos;
        std::vector<double> z_pos;
        std::vector<double> y_radius;
        std::vector<double> r_values;
    };
}

static WheelSetup setupWheel(const Profile &wp, const double r0, const double b0,
                             const double state_y, const double state_z, const Side side) {
    WheelSetup ws;

    // Reference index: first point where y >= 0 (tread center)
    size_t ref_idx = 0;
    for (size_t i = 0; i < wp.y.size(); ++i) {
        if (wp.y[i] >= 0.0) {
            ref_idx = i;
            break;
        }
    }
    const double z_ref = wp.z[ref_idx];

    const size_t n = wp.y.size();
    ws.y_pos.resize(n);
    ws.z_pos.resize(n);
    ws.y_radius.resize(n);
    ws.r_values.resize(n);

    if (side == Side::RIGHT) {
        for (size_t i = 0; i < n; ++i) {
            ws.y_pos[i] = wp.y[i] + b0 + state_y;
            ws.z_pos[i] = wp.z[i] + state_z;
            ws.y_radius[i] = wp.y[i] + b0 + state_y;
            ws.r_values[i] = r0 + (wp.z[i] - z_ref);
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            size_t j = n - 1 - i;
            ws.y_pos[i] = -wp.y[j] - b0 + state_y;
            ws.z_pos[i] = wp.z[j] + state_z;
            ws.y_radius[i] = -wp.y[j] - b0 + state_y;
            ws.r_values[i] = r0 + (wp.z[j] - z_ref);
        }
    }

    return ws;
}

namespace {
    // ---------------------------------------------------------
    // Rail setup — matches Python rail class
    // ---------------------------------------------------------
    struct RailSetup {
        std::vector<double> y_pos;
        std::vector<double> z_pos;
    };
}

static RailSetup setupRail(const Profile &rp, double gauge, int inclination, Side side) {
    RailSetup rs;

    const double rail_cant = -std::atan(1.0 / inclination);

    std::vector<double> y_rot(rp.y.size());
    std::vector<double> z_rot(rp.z.size());

    for (size_t i = 0; i < rp.y.size(); ++i) {
        y_rot[i] = rp.y[i];
        z_rot[i] = rp.z[i];
        rotatePoint(y_rot[i], z_rot[i], rail_cant);
    }

    // Gauge corner: first index where z_rotated <= 0.014
    size_t gc_idx = 0;
    for (size_t i = 0; i < z_rot.size(); ++i) {
        if (z_rot[i] <= 0.014) {
            gc_idx = i;
            break;
        }
    }

    const size_t n = y_rot.size();
    rs.y_pos.resize(n);
    rs.z_pos.resize(n);

    if (side == Side::RIGHT) {
        for (size_t i = 0; i < n; ++i) {
            rs.y_pos[i] = y_rot[i] - y_rot[gc_idx] + gauge / 2.0;
            rs.z_pos[i] = z_rot[i];
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            size_t j = n - 1 - i;
            rs.y_pos[i] = -y_rot[j] + y_rot[gc_idx] - gauge / 2.0;
            rs.z_pos[i] = z_rot[j];
        }
    }

    return rs;
}

namespace {
    // ---------------------------------------------------------
    // Interpolate wheel z and radius onto rail y-coordinates
    // ---------------------------------------------------------
    struct InterpData {
        std::vector<double> wheel_z;
        std::vector<double> radius;
    };
}

static InterpData interpolateWheelOnRail(const WheelSetup &ws, const RailSetup &rs) {
    InterpData data;
    data.wheel_z.resize(rs.y_pos.size());
    data.radius.resize(rs.y_pos.size());

    for (size_t i = 0; i < rs.y_pos.size(); ++i) {
        data.wheel_z[i] = linearInterp(ws.y_pos, ws.z_pos, rs.y_pos[i]);
        data.radius[i] = linearInterp(ws.y_radius, ws.r_values, rs.y_pos[i]);
    }
    return data;
}

// ---------------------------------------------------------
// Find contact regions: where wheel_z - rail_z > 0
// ---------------------------------------------------------
static std::vector<std::pair<size_t, size_t> > findContactRegions(const std::vector<double> &wheel_z,
                                                                  const std::vector<double> &rail_z) {
    std::vector<size_t> contact_idx;
    for (size_t i = 0; i < wheel_z.size() && i < rail_z.size(); ++i) {
        if (wheel_z[i] - rail_z[i] > constants::CONTACT_TOL)
            contact_idx.push_back(i);
    }
    if (contact_idx.empty()) return {};

    std::vector<std::pair<size_t, size_t> > regions;
    size_t start = contact_idx[0];
    size_t end = contact_idx[0];

    for (size_t i = 1; i < contact_idx.size(); ++i) {
        if (contact_idx[i] == end + 1) {
            end = contact_idx[i];
        } else {
            regions.emplace_back(start, end);
            start = end = contact_idx[i];
        }
    }
    regions.emplace_back(start, end);

    std::ranges::sort(regions, [](const auto &a, const auto &b) {
        return (b.second - b.first) < (a.second - a.first);
    });

    return regions;
}

namespace {
    // ---------------------------------------------------------
    // Process a contact region: rotate to normal direction,
    // re-discretize to uniform grid (Eq. 4-7)
    // ---------------------------------------------------------
    struct ContactPatch {
        std::vector<double> local_x;
        std::vector<double> wheel_z_norm;
        std::vector<double> rail_z_norm;
        std::vector<double> radius_norm;
        double rotation_angle{};
    };
}

// ---------------------------------------------------------
// Process a contact region: optimized allocation footprint
// ---------------------------------------------------------
static ContactPatch processContactRegion(
    const RailSetup &rs,
    const InterpData &id,
    size_t start_idx,
    size_t end_idx,
    int discretization) {
    ContactPatch cp;
    const size_t len = end_idx - start_idx + 1;

    // Use thread_local static transient arrays to entirely eliminate inner-loop allocations
    thread_local std::vector<double> rail_y_c;
    thread_local std::vector<double> rail_z_c;
    thread_local std::vector<double> wheel_z_c;
    thread_local std::vector<double> radius_c;
    thread_local std::vector<double> y_t, z_rail_t, z_wheel_t, z_rad_t;
    thread_local std::vector<double> y_rot, z_rail_rot, z_wheel_rot, z_rad_rot;

    rail_y_c.resize(len);
    rail_z_c.resize(len);
    wheel_z_c.resize(len);
    radius_c.resize(len);
    y_t.resize(len);
    z_rail_t.resize(len);
    z_wheel_t.resize(len);
    z_rad_t.resize(len);
    y_rot.resize(len);
    z_rail_rot.resize(len);
    z_wheel_rot.resize(len);
    z_rad_rot.resize(len);

    for (size_t i = 0; i < len; ++i) {
        size_t idx = start_idx + i;
        rail_y_c[i] = rs.y_pos[idx];
        rail_z_c[i] = rs.z_pos[idx];
        wheel_z_c[i] = id.wheel_z[idx];
        radius_c[i] = id.radius[idx];
    }

    // Eq. 4: contact rotation angle alpha
    double dy = rail_y_c.back() - rail_y_c.front();
    double dz = rail_z_c.back() - rail_z_c.front();
    double alpha = -std::atan2(dz, dy);
    cp.rotation_angle = alpha;

    double ca = std::cos(alpha);
    double sa = std::sin(alpha);

    // Translate to origin (Eq. 6 normalization)
    for (size_t i = 0; i < len; ++i) {
        y_t[i] = rail_y_c[i] - rail_y_c.front();
        z_rail_t[i] = rail_z_c[i] - rail_z_c.front();
        z_wheel_t[i] = wheel_z_c[i] - wheel_z_c.front();
        z_rad_t[i] = radius_c[i] - radius_c.front();
    }

    // Apply rotation R_alpha (Eq. 6)
    for (size_t i = 0; i < len; ++i) {
        y_rot[i] = ca * y_t[i] - sa * z_rail_t[i];
        z_rail_rot[i] = sa * y_t[i] + ca * z_rail_t[i];
        z_wheel_rot[i] = sa * y_t[i] + ca * z_wheel_t[i];
        z_rad_rot[i] = sa * y_t[i] + ca * z_rad_t[i];
    }

    // Re-discretize onto uniform grid (Eq. 7)
    cp.local_x.resize(static_cast<size_t>(discretization));
    double y0 = y_rot.front();
    double y1 = y_rot.back();
    for (int i = 0; i < discretization; ++i) {
        cp.local_x[static_cast<size_t>(i)] = y0 + static_cast<double>(i) * (y1 - y0)
                                             / static_cast<double>(discretization - 1);
    }

    // Interpolate onto uniform grid using PCHIP
    geometry::PCHIPInterpolator rail_interp(y_rot, z_rail_rot);
    geometry::PCHIPInterpolator wheel_interp(y_rot, z_wheel_rot);
    geometry::PCHIPInterpolator radius_interp(y_rot, z_rad_rot);
    cp.rail_z_norm = rail_interp.evaluateBatch(cp.local_x);
    cp.wheel_z_norm = wheel_interp.evaluateBatch(cp.local_x);
    cp.radius_norm = radius_interp.evaluateBatch(cp.local_x);

    // Restore vertical offset (Eq. 7): add back initial values
    for (size_t i = 0; i < cp.wheel_z_norm.size(); ++i) {
        cp.wheel_z_norm[i] += wheel_z_c.front();
        cp.rail_z_norm[i] += rail_z_c.front();
        cp.radius_norm[i] += radius_c.front();
    }
    return cp;
}

namespace {
    // ---------------------------------------------------------
    // Contact result for single lateral position
    // ---------------------------------------------------------
    struct ContactResultSingle {
        double state_y; // lateral displacement [m]
        double state_z; // solved vertical position [m]
        double normal_force; // [N]
        double q_force; // [N]
        double y_force; // [N]
        double approach; // equivalent penetration [m]
        double semi_axis_a; // [m]
        double semi_axis_b; // [m]
        double a_b_ratio; // [-]
        double contact_angle; // [rad]
        double centroid_y; // [m]
        double contact_pos_rail; // [mm] relative to rail center
        double contact_pos_wheel; // [mm] relative to wheel center
        double patch_length; // [mm] — 2*a
        double patch_width; // [mm] — 2*b
        double patch_area; // [mm²]
        bool contact_found;
    };
}

// ---------------------------------------------------------------
// Compute contact for given state_y, state_z
// Encapsulates the full pipeline: wheel setup → interpolation →
// contact detection → rotation → solving → position calc (Eq. 23)
// ---------------------------------------------------------------
static std::optional<ContactResultSingle> computeContactOnce(
    const Profile &wheel_prof, const RailSetup &rs, double r0,
    double b0,
    double gauge, double state_y, double state_z, Side side,
    const MaterialProperties &material, int discretization) {
    WheelSetup ws = setupWheel(wheel_prof, r0, b0, state_y, state_z, side);
    InterpData id = interpolateWheelOnRail(ws, rs);

    auto regions = findContactRegions(id.wheel_z, rs.z_pos);

    if (regions.empty()) return std::nullopt; // TRUE NO CONTACT (Wheel is floating)

    // Process up to 3 patches (matches Python n_patches=3)
    constexpr int MAX_PATCHES = 3;
    int n_patches = static_cast<int>(std::min(regions.size(), static_cast<size_t>(MAX_PATCHES)));

    double total_q = 0.0;
    double total_normal = 0.0;
    double total_y = 0.0;
    double total_patch_area = 0.0; // Cumulative area tracker

    // Store dominant patch (largest region) for geometry reporting
    std::optional<ContactResultSingle> dominant;
    double max_patch_area = -1.0; // Track the largest area found

    for (int p = 0; p < n_patches; ++p) {
        auto [start_idx, end_idx] = regions[static_cast<size_t>(p)];

        // Initial values for contact position back-transform (Eq. 23)
        double y_rail_ini = rs.y_pos[start_idx];
        double z_rail_ini = rs.z_pos[start_idx];
        double z_wheel_ini = id.wheel_z[start_idx];

        auto cp = processContactRegion(rs, id, start_idx, end_idx, discretization);

        // Eq. 14: local rolling radius at contact center
        size_t mid = cp.radius_norm.size() / 2;
        double r_local = cp.radius_norm[mid] / std::cos(cp.rotation_angle);

        ContactResult result = contact::EquivalentContactSolver::solveContact(
            cp.local_x, cp.wheel_z_norm, cp.rail_z_norm, r_local, material, cp.rotation_angle);

        if (!result.contact_found) continue;

        // 1. Accumulate total physical force components across all active patches
        total_q += result.q_force;
        total_normal += result.normal_force;
        total_y += result.y_force;

        // 2. Calculate the individual patch area and add to the total aggregate
        double current_patch_area_m2 = std::numbers::pi * result.semi_axis_a * result.semi_axis_b;
        total_patch_area += current_patch_area_m2 * 1e6;

        // 3. Track geometric dominance based on individual patch size
        if (current_patch_area_m2 > max_patch_area) {
            max_patch_area = current_patch_area_m2;

            // Contact position back-transform (Eq. 23)
            double cy = result.centroid_y;
            double alpha = cp.rotation_angle;
            double ca = std::cos(alpha);
            double sa = std::sin(alpha);

            double z_rail_at_cy = linearInterp(cp.local_x, cp.rail_z_norm, cy);
            double z_wheel_at_cy = linearInterp(cp.local_x, cp.wheel_z_norm, cy);

            // Eq. 23a: rail contact position
            double y_rail_cy = ca * cy + sa * (z_rail_at_cy - z_rail_ini) + y_rail_ini;
            // Eq. 23b: wheel contact position
            double y_wheel_cy = ca * cy + sa * (z_wheel_at_cy - z_wheel_ini) + y_rail_ini;

            // Relative to rail/wheel center
            double contact_pos_rail, contact_pos_wheel;
            if (side == Side::LEFT) {
                contact_pos_rail = (y_rail_cy + gauge / 2.0) * 1000.0;
                contact_pos_wheel = (y_wheel_cy + b0) * 1000.0;
            } else {
                contact_pos_rail = (y_rail_cy - gauge / 2.0) * 1000.0;
                contact_pos_wheel = (y_wheel_cy - b0) * 1000.0;
            }

            ContactResultSingle sr{};
            sr.state_y = state_y;
            sr.state_z = state_z;
            sr.normal_force = 0; // filled below with totals
            sr.q_force = 0;
            sr.y_force = 0;
            sr.approach = result.approach;
            sr.semi_axis_a = result.semi_axis_a;
            sr.semi_axis_b = result.semi_axis_b;
            sr.a_b_ratio = result.a_b_ratio;
            sr.contact_angle = result.contact_angle;
            sr.centroid_y = result.centroid_y;
            sr.contact_pos_rail = contact_pos_rail;
            sr.contact_pos_wheel = contact_pos_wheel;
            sr.patch_length = result.semi_axis_a * 2000.0;
            sr.patch_width = result.semi_axis_b * 2000.0;
            // sr.patch_area = current_patch_area * 1e6; // Convert to mm²
            sr.patch_area = 0.0; // Overwritten below with cumulative sum
            sr.contact_found = true;
            dominant = sr; // Captured the true dominant geometrical properties
        }
    }

    if (!dominant) {
        // High-penetration solver failure fallback:
        // Profiles intersect (regions exist) but inner solver couldn't handle the size.
        ContactResultSingle failure_fallback{};
        failure_fallback.state_y = state_y;
        failure_fallback.state_z = state_z;
        failure_fallback.q_force = constants::SOLVER_FAILURE_FORCE; // Assign a massive dummy force
        failure_fallback.normal_force = constants::SOLVER_FAILURE_FORCE;
        failure_fallback.contact_found = false;
        return failure_fallback;
    }

    // Replace single-patch forces with SUMMED forces
    dominant->normal_force = total_normal;
    dominant->q_force = total_q;
    dominant->y_force = total_y;
    dominant->patch_area = total_patch_area; // <-- OVERWRITE WITH SUMMED PATCH AREAS FOR PLOTTING

    return dominant;
}

// =========================================================
// Core Logic: Runs the entire simulation sweep
// =========================================================
static int runContactSweep(
    const std::string &wheelProfilePath,
    const std::string &railProfilePath,
    double b0,
    double r0,
    double gauge,
    int railInclination,
    double targetQ,
    Side side,
    double dyMin,
    double dyMax,
    double dyStep,
    double qTol,
    int maxIter,
    double zLoInit,
    double zHiInit,
    double dzInitialGuess,
    const MaterialProperties &material) {
    try {
        Profile raw_wheel = geometry::loadProfileFromFile(wheelProfilePath);
        Profile raw_rail = geometry::loadProfileFromFile(railProfilePath);
        Profile wheel_prof = resampleProfile(raw_wheel);
        Profile rail_prof = resampleProfile(raw_rail);
        std::cout << "Loaded profiles:\n";
        std::cout << "  Wheel: " << wheel_prof.y.size() << " points\n";
        std::cout << "  Rail:  " << rail_prof.y.size() << " points\n\n";
        RailSetup rs = setupRail(rail_prof, gauge, railInclination, side);
        std::cout << "Sweep configuration:\n";
        std::cout << "  Lateral range:       [" << dyMin * 1000.0 << ", " << dyMax * 1000.0 << "] mm\n";
        std::cout << "  Lateral step:        " << dyStep * 1000.0 << " mm\n";
        std::cout << "  Target Q force:      " << targetQ / 1000.0 << " kN\n";
        std::cout << "  Side:                " << (side == Side::LEFT ? "LEFT" : "RIGHT") << "\n\n";

#ifdef USE_OPENMP
        const int num_threads = omp_get_max_threads();
        std::cout << " OpenMP threads: " << num_threads << "\n";
#else
        std::cout << " OpenMP: Disabled (single-threaded)\n";
#endif
        std::cout << "\n";

        // --- Set CSV file ---
        const std::string outputPath = "sweep_results.csv";
        std::ofstream csv(outputPath);
        csv << std::fixed << std::setprecision(6);
        csv << "dy_mm,state_z_mm,normal_force_N,q_force_N,y_force_N,"
                << "semi_axis_a_mm,semi_axis_b_mm,ab_ratio,"
                << "patch_length_mm,patch_width_mm,patch_area_mm2,"
                << "approach_mm,contact_angle_deg,centroid_y_mm,"
                << "contact_pos_rail_mm,contact_pos_wheel_mm\n";

        auto t0_total = std::chrono::high_resolution_clock::now();

        // --- Sweeping over lateral positions ---
        // ---------------------------------------
        // 1. Calculate total steps: (Span / Step) + 1 for an inclusive (<=) loop
        // std::round handles any tiny floating-point inaccuracies in the division

        // --- Pre-allocate results vector (thread-safe indexing) ---
        const size_t total_steps = static_cast<size_t>(std::round((dyMax - dyMin) / dyStep)) + 1;
        std::vector<std::optional<ContactResultSingle> > results(total_steps);
        double dz_warm = dzInitialGuess; // Warm-start for next position

#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic, 8) // NOLINT(*-use-default-none)
        // 2. Using size_t induction variable to satisfy Clang-Tidy
        for (size_t i = 0; i < total_steps; ++i) {
#else
            // 2. Use a size_t induction variable to satisfy Clang-Tidy
            for (size_t i = 0; i < total_steps; ++i) {
#endif
            // 3. Compute state_y fresh each iteration to prevent rounding drift
            double state_y = dyMin + (static_cast<double>(i) * dyStep);

            // Expand limits to allow negative values (lifting the wheel)
            double z_lo = zLoInit;
            double z_hi = std::max(dz_warm * 3.0, zHiInit);

            auto t0 = std::chrono::high_resolution_clock::now();

            std::optional<ContactResultSingle> best;
            double best_q = 0.0;

            // Bisection loop
            for (int iter = 0; iter < maxIter; ++iter) {
                double z_mid = (z_lo + z_hi) / 2.0;

                auto result = computeContactOnce(
                    wheel_prof, rs, r0, b0, gauge,
                    state_y, z_mid, side, material,
                    constants::DEFAULT_DISCRETIZATION);

                if (!result) {
                    // No contact entirely -> wheel is floating too high up, push it down
                    z_lo = z_mid;
                    continue;
                }

                best = result;
                best_q = result->q_force;

                if (std::abs(best_q - targetQ) < qTol) {
                    break;
                }

                if (best_q < targetQ) {
                    z_lo = z_mid; // Force too low -> push wheel deeper down
                } else {
                    z_hi = z_mid; // Force too high -> lift wheel higher up
                }
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            auto us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);

            // Store result at thread-local index (no lock needed)
            results[i] = best;
            (void) us; // Suppress unused warning in release builds
        }

        // --- Sequential CSV write after all threads join ---
        int count_converged = 0;
        for (size_t i = 0; i < total_steps; ++i) {
            double state_y = dyMin + static_cast<double>(i) * dyStep;

            if (results[i] && std::abs(results[i]->q_force - targetQ) < qTol * 5) {
                csv << state_y * 1000.0 << ","
                        << results[i]->state_z * 1000.0 << ","
                        << results[i]->normal_force << ","
                        << results[i]->q_force << ","
                        << results[i]->y_force << ","
                        << results[i]->semi_axis_a * 1000.0 << ","
                        << results[i]->semi_axis_b * 1000.0 << ","
                        << results[i]->a_b_ratio << ","
                        << results[i]->patch_length << ","
                        << results[i]->patch_width << ","
                        << results[i]->patch_area << ","
                        << results[i]->approach * 1000.0 << ","
                        << results[i]->contact_angle * 180.0 / M_PI << ","
                        << results[i]->centroid_y * 1000.0 << ","
                        << results[i]->contact_pos_rail << ","
                        << results[i]->contact_pos_wheel << "\n";
                ++count_converged;
            } else {
                // No contact or no convergence — write zeros
                csv << state_y * 1000.0 << ",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0\n";
            }
        }

        csv.close();

        // --- Progress report ---
        auto t1_total = std::chrono::high_resolution_clock::now();
        auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1_total - t0_total);

        std::cout << "=== Sweep Complete ===\n";
        std::cout << "Total positions: " << total_steps << "\n";
        std::cout << "Converged: " << count_converged << " ("
                << std::fixed << std::setprecision(1)
                << (100.0 * count_converged / static_cast<int>(total_steps)) << "%)\n";
        std::cout << "Elapsed: " << total_ms.count() << " ms\n";

#ifdef USE_OPENMP
        std::cout << "Threads: " << omp_get_max_threads() << "\n";
        std::cout << "Avg per step: " << (total_ms.count() / std::max(1, count_converged))
                << " ms (parallel)\n";
#else
        std::cout << "Avg per step: " << (total_ms.count() / std::max(1, count_converged))
                << " ms (serial)\n";
#endif

        std::cout << "Results written to " << outputPath << "\n";

        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

// =========================================================
// Entry point: Inputs and Output only
// =========================================================
int main() {
    std::cout << "=== Equivalent Elastic Contact (Lateral Sweep) ===\n\n";

    // --- Input Configuration ---
    const std::string wheelProfilePath = "S1002.txt";
    const std::string railProfilePath = "UIC60.txt";

    constexpr double b0 = 0.75; // Semi-axis distance [m]
    constexpr double r0 = 0.46; // Nominal radius [m]
    constexpr double gauge = 1.435; // Track gauge [m]
    constexpr int rail_inclination = 40; // 1:40
    constexpr double target_q = 50000.0; // Target Q force [N]
    constexpr auto side = Side::RIGHT;

    constexpr double dy_min = -19.5e-3; // [mm]
    constexpr double dy_max = 19.5e-3; // [mm]
    constexpr double dy_step = 0.1e-3; // [mm]

    constexpr double q_tol = 50.0; // Q force tolerance [N]
    constexpr int max_iter = 60; // Bisection iterations

    constexpr double z_lo_init = -0.5; // Allows the wheel to climb up (flange zone) [m]
    constexpr double z_hi_init = 0.05; // Upper bound for state_z [m]
    constexpr double dz_initial_guess = 0.002;

    MaterialProperties material;
    material.elastic_modulus = 210e9; // [N/m2]
    material.shear_coefficient = 0.28;

    // Execute the simulation logic
    return runContactSweep(
        wheelProfilePath, railProfilePath,
        b0, r0, gauge, rail_inclination, target_q,
        side, dy_min, dy_max, dy_step,
        q_tol, max_iter, z_lo_init, z_hi_init,
        dz_initial_guess, material);
}
