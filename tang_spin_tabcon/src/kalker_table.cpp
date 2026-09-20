#include "kalker_table.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

namespace {
    /**
 * @struct InternalEntry
 * @brief Internal representation of one stored table point (normalized forces).
 */
    struct InternalEntry {
        double Fx, Fy, Mz;
    };
}

// ----- con93 tabcon.dat true layout structure (itab = 7) -----
// Both a>b and +/-phi are stored explicitly, so lookups need no axis swap;
// Only upsx/upsy (stored non-negative) require sign reduction.
namespace {
    constexpr size_t ITAB = 7;
    constexpr size_t NAOB = 2 * ITAB;
    constexpr size_t NUPS = 2 * (ITAB + 1);
    constexpr size_t NPHI = 2 * NUPS;
    constexpr size_t NSUB = NUPS * NUPS;

    // ups grid values (upsx and upsy share 16-point grid).
    std::vector<double> make_ups() {
        std::vector<double> v(NUPS);
        for (size_t i = 0; i <= ITAB; ++i) v[i] = static_cast<double>(i) / static_cast<double>(ITAB); // 0..1
        for (size_t i = 0; i <= ITAB; ++i)
            v[ITAB + 1 + i] = static_cast<double>(ITAB) / std::max(
                                  1e-5, static_cast<double>(i)); // inf..1
        return v;
    }

    // phi grid values: ups list, then its negation.
    std::vector<double> make_phi() {
        const std::vector<double> u = make_ups();
        std::vector<double> v(NPHI);
        for (size_t i = 0; i < NUPS; ++i) {
            v[i] = u[i];
            v[NUPS + i] = -u[i];
        }
        return v;
    }

    // aob grid values in file-block order.
    std::vector<double> make_aob() {
        std::vector<double> v(NAOB);
        for (size_t i = 1; i <= ITAB; ++i) v[i - 1] = static_cast<double>(i) / static_cast<double>(ITAB);
        // 1/7..1 (a<=b)
        for (size_t i = 1; i <= ITAB; ++i) v[ITAB + i - 1] = static_cast<double>(ITAB) / static_cast<double>(i);
        // 7..1  (a>b)
        return v;
    }

    /**
     * @struct Axis
     * @brief A lookup axis whose stored sample points are not necessarily sorted
     *        (the ups/phi/aob grids fold a "linear" and a "reciprocal" branch).
     *
     * Holds the sample values sorted ascending together with the data index each
     * value maps to, so interpolation can bracket any query monotonically.
     */
    struct Axis {
        std::vector<double> val; // sorted ascending
        std::vector<size_t> idx; // val[k] lives at data column idx[k]

        void build(const std::vector<double> &raw) {
            std::vector<std::pair<double, size_t> > p;
            p.reserve(raw.size());
            for (size_t i = 0; i < raw.size(); ++i) p.emplace_back(raw[i], i);
            std::ranges::sort(p, [](const auto &a, const auto &b) { return a.first < b.first; });
            // Drop duplicate values (e.g. aob=1 and ups=1 appear twice with identical data).
            for (const auto &[v, i]: p) {
                if (!val.empty() && std::abs(v - val.back()) < 1e-12) continue;
                val.push_back(v);
                idx.push_back(i);
            }
        }

        // Bracket a (clamped) query: returns lower sorted-slot k and weight w
        // so the result is corner[idx[k]]*(1-w) + corner[idx[k+1]]*w.
        void locate(double q, size_t &k, double &w) const {
            q = std::clamp(q, val.front(), val.back());
            const auto it = std::ranges::upper_bound(val, q);
            auto hi = static_cast<size_t>(std::distance(val.begin(), it));
            if (hi == 0) hi = 1;
            if (hi >= val.size()) hi = val.size() - 1;
            k = hi - 1;
            const double lo_v = val[k], hi_v = val[hi];
            w = (hi_v - lo_v < 1e-15) ? 0.0 : (q - lo_v) / (hi_v - lo_v);
        }
    };
}

/**
 * @struct KalkerTable::Impl
 */
struct KalkerTable::Impl {
    std::vector<InternalEntry> data;
    Axis ax_aob, ax_ups, ax_phi;

    Impl() {
        ax_aob.build(make_aob());
        ax_ups.build(make_ups());
        ax_phi.build(make_phi());
    }

    [[nodiscard]] const InternalEntry &at(const size_t b, const size_t iupsx, const size_t iupsy,
                                          const size_t iphi) const {
        return data[((b * NUPS + iupsx) * NUPS + iupsy) * NPHI + iphi];
    }
};

KalkerTable::KalkerTable() : pimpl(std::make_unique<Impl>()) {
}

KalkerTable::~KalkerTable() = default;

bool KalkerTable::load(const std::string &filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    for (int i = 0; i < 7; ++i) {
        if (!std::getline(file, line)) return false; // skip 7-line header
    }

    pimpl->data.assign(NAOB * NSUB * NPHI, InternalEntry{});

    // Read the sub-blocks strictly in file order; the first two integers
    // on each "a/b" header line are positional bookkeeping, ignore them.
    size_t written = 0;
    const size_t total = pimpl->data.size();
    while (written < total && std::getline(file, line)) {
        if (line.find("a/b") == std::string::npos) continue;
        for (size_t k = 0; k < NPHI; ++k) {
            if (!(file >> pimpl->data[written].Fx
                  >> pimpl->data[written].Fy
                  >> pimpl->data[written].Mz)) {
                return false;
            }
            ++written;
        }
        std::getline(file, line); // consume the rest of the last entry line
    }

    return written == total;
}

KalkerResult KalkerTable::lookup(const double aob, const double ux, const double uy, const double phi) const {
    // Reduce to the stored octant: upsx, upsy >= 0. The mirror relations
    // (reflection about the x-z and y-z planes) map the sign of ux and uy
    // onto the output forces and flip phi (crucially):
    // Fx = sx * Fx_s ,  Fy = sy * Fy_s ,  Mz = sx*sy * Mz_s
    // evaluated at ( |ux|, |uy|, phi_eff = sx*sy*phi ).
    const double sx = (ux < 0) ? -1.0 : 1.0;
    const double sy = (uy < 0) ? -1.0 : 1.0;
    const double axu = std::abs(ux);
    const double ayu = std::abs(uy);
    const double phi_eff = sx * sy * phi;

    size_t ba = 0, bx = 0, by = 0, bp = 0;
    double wa = 0.0, wx = 0.0, wy = 0.0, wp = 0.0;
    pimpl->ax_aob.locate(aob, ba, wa);
    pimpl->ax_ups.locate(axu, bx, wx);
    pimpl->ax_ups.locate(ayu, by, wy);
    pimpl->ax_phi.locate(phi_eff, bp, wp);

    // aob / ups / phi sorted-slot -> data indices for the two bracketing samples.
    const size_t ai[2] = {pimpl->ax_aob.idx[ba], pimpl->ax_aob.idx[ba + 1]};
    const size_t xi[2] = {pimpl->ax_ups.idx[bx], pimpl->ax_ups.idx[bx + 1]};
    const size_t yi[2] = {pimpl->ax_ups.idx[by], pimpl->ax_ups.idx[by + 1]};
    const size_t pi[2] = {pimpl->ax_phi.idx[bp], pimpl->ax_phi.idx[bp + 1]};
    const double aw[2] = {1.0 - wa, wa};
    const double xw[2] = {1.0 - wx, wx};
    const double yw[2] = {1.0 - wy, wy};
    const double pw[2] = {1.0 - wp, wp};

    // Quadrilinear blend over the 16 surrounding corners.
    InternalEntry acc{0.0, 0.0, 0.0};
    for (size_t a = 0; a < 2; ++a)
        for (size_t x = 0; x < 2; ++x)
            for (size_t y = 0; y < 2; ++y)
                for (size_t p = 0; p < 2; ++p) {
                    const double w = aw[a] * xw[x] * yw[y] * pw[p];
                    const InternalEntry &e = pimpl->at(ai[a], xi[x], yi[y], pi[p]);
                    acc.Fx += w * e.Fx;
                    acc.Fy += w * e.Fy;
                    acc.Mz += w * e.Mz;
                }

    KalkerResult res{};
    res.Fx = sx * acc.Fx;
    res.Fy = sy * acc.Fy;
    res.Mz = sx * sy * acc.Mz;
    return res;
}
