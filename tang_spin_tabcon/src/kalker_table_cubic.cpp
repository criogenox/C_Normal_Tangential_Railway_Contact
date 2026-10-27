#include "kalker_table_cubic.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

// This library shares tabcon.dat's layout, grids, and symmetry rules with the
// linear KalkerTable. It is intentionally self-contained (loader + grid construction
// duplicated) so it depends on nothing from the linear library. The ONLY conceptual
// difference is the interpolation kernel: monotone cubic Hermite (PCHIP/Fritsch-Carlson)
// instead of multilinear.

/**
 * @struct InternalEntry
 * @brief One stored table point (normalized forces).
 */
struct InternalEntry {
    double Fx, Fy, Mz;
};

namespace {
    constexpr size_t ITAB = 7;
    constexpr size_t NAOB = 2 * ITAB;
    constexpr size_t NUPS = 2 * (ITAB + 1);
    constexpr size_t NPHI = 2 * NUPS;
    constexpr size_t NSUB = NUPS * NUPS;

    // ups grid values (upsx and upsy share this 16-point grid).
    std::vector<double> make_ups() {
        std::vector<double> v(NUPS);
        for (size_t i = 0; i <= ITAB; ++i) v[i] = static_cast<double>(i) / ITAB;
        for (size_t i = 0; i <= ITAB; ++i) v[ITAB + 1 + i] = ITAB / std::max(1e-5, static_cast<double>(i));
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
        for (size_t i = 1; i <= ITAB; ++i) v[ITAB + i - 1] = static_cast<double>(ITAB) / static_cast<double>(i);
        return v;
    }

    /**
     * @struct Axis
     * @brief Grid samples sorted ascending, mapped back to their data column, so
     *        any query can be bracketed monotonically over the folded grid.
     */
    struct Axis {
        std::vector<double> val; // sorted ascending, de-duplicated
        std::vector<size_t> idx; // val[k] lives at data column idx[k]

        void build(const std::vector<double> &raw) {
            std::vector<std::pair<double, size_t> > p;
            p.reserve(raw.size());
            for (size_t i = 0; i < raw.size(); ++i) p.emplace_back(raw[i], i);
            std::ranges::sort(p, [](const auto &a, const auto &b) { return a.first < b.first; });
            for (const auto &[v, i]: p) {
                if (!val.empty() && std::abs(v - val.back()) < 1e-12) continue;
                val.push_back(v);
                idx.push_back(i);
            }
        }

        // Lower node index k of the interval [val[k], val[k+1]] containing (clamped) q.
        [[nodiscard]] int interval(double q) const {
            q = std::clamp(q, val.front(), val.back());
            const auto it = std::ranges::upper_bound(val, q);
            int hi = static_cast<int>(std::distance(val.begin(), it));
            if (hi <= 0) hi = 1;
            if (hi >= static_cast<int>(val.size())) hi = static_cast<int>(val.size()) - 1;
            return hi - 1;
        }

        [[nodiscard]] int size() const { return static_cast<int>(val.size()); }
    };

    /**
     * @struct Stencil
     * @brief Up to 4 consecutive axis nodes around an interval, with their axis
     *        coordinates (x) and data columns (col). Clamped and de-duplicated at
     *        the table boundaries, so n is 2..4.
     */
    struct Stencil {
        size_t n = 0;
        double x[4]{};
        size_t col[4]{};
    };

    Stencil make_stencil(const Axis &ax, const int k) {
        const int raw[4] = {k - 1, k, k + 1, k + 2};
        const int sz = ax.size();
        Stencil s;
        int last = -1;
        for (const int c: raw) {
            const int pos = std::clamp(c, 0, sz - 1);
            if (pos == last) continue; // drop clamp-induced repeats at the ends
            s.x[s.n] = ax.val[static_cast<size_t>(pos)];
            s.col[s.n] = ax.idx[static_cast<size_t>(pos)];
            ++s.n;
            last = pos;
        }
        return s;
    }

    // Fritsch-Carlson weighted-harmonic-mean slope at an interior node.
    // dl,dr: secant slopes of the left/right intervals; hl,hr: their widths.
    // Returns 0 at local extrema (opposite-sign or zero secants) to preserve monotonicity.
    double pchip_slope(const double dl, const double dr, const double hl, const double hr) {
        if (dl * dr <= 0.0) return 0.0;
        const double w1 = 2.0 * hr + hl;
        const double w2 = hr + 2.0 * hl;
        return (w1 + w2) / (w1 / dl + w2 / dr);
    }

    // Evaluate a monotone cubic Hermite through the (x[i], y[i]) stencil at q.
    // Slopes: interior nodes use pchip_slope; boundary nodes use the one-sided
    // secant (safe, non-overshooting). Reduces to exact linear when m == 2.
    // Because a 4-node local stencil around interval k carries exactly the nodes
    // PCHIP needs for the slopes at k and k+1, this matches full-axis PCHIP there.
    InternalEntry mono_cubic(const double *x, const InternalEntry *y, const size_t m, const double q) {
        if (m == 1) return y[0];

        double h[3];
        InternalEntry d[3];
        for (size_t i = 0; i < m - 1; ++i) {
            h[i] = x[i + 1] - x[i];
            const double inv = (h[i] > 1e-30) ? 1.0 / h[i] : 0.0;
            d[i].Fx = (y[i + 1].Fx - y[i].Fx) * inv;
            d[i].Fy = (y[i + 1].Fy - y[i].Fy) * inv;
            d[i].Mz = (y[i + 1].Mz - y[i].Mz) * inv;
        }

        InternalEntry mslope[4];
        for (size_t k = 0; k < m; ++k) {
            if (k == 0) {
                mslope[0] = d[0];
            } else if (k == m - 1) {
                mslope[k] = d[m - 2];
            } else {
                mslope[k].Fx = pchip_slope(d[k - 1].Fx, d[k].Fx, h[k - 1], h[k]);
                mslope[k].Fy = pchip_slope(d[k - 1].Fy, d[k].Fy, h[k - 1], h[k]);
                mslope[k].Mz = pchip_slope(d[k - 1].Mz, d[k].Mz, h[k - 1], h[k]);
            }
        }

        size_t j = 0;
        while (j < m - 2 && q > x[j + 1]) ++j;
        const double hh = x[j + 1] - x[j];
        const double t = (hh > 1e-30) ? (q - x[j]) / hh : 0.0;
        const double t2 = t * t, t3 = t2 * t;
        const double h00 = 2 * t3 - 3 * t2 + 1;
        const double h10 = t3 - 2 * t2 + t;
        const double h01 = -2 * t3 + 3 * t2;
        const double h11 = t3 - t2;

        InternalEntry r{};
        r.Fx = h00 * y[j].Fx + h10 * hh * mslope[j].Fx + h01 * y[j + 1].Fx + h11 * hh * mslope[j + 1].Fx;
        r.Fy = h00 * y[j].Fy + h10 * hh * mslope[j].Fy + h01 * y[j + 1].Fy + h11 * hh * mslope[j + 1].Fy;
        r.Mz = h00 * y[j].Mz + h10 * hh * mslope[j].Mz + h01 * y[j + 1].Mz + h11 * hh * mslope[j + 1].Mz;
        return r;
    }
}

/**
 * @struct KalkerTableCubic::Impl
 */
struct KalkerTableCubic::Impl {
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

KalkerTableCubic::KalkerTableCubic() : pimpl(std::make_unique<Impl>()) {
}

KalkerTableCubic::~KalkerTableCubic() = default;

bool KalkerTableCubic::load(const std::string &filename) const {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    for (int i = 0; i < 7; ++i) {
        if (!std::getline(file, line)) return false; // skip 7-line header
    }

    pimpl->data.assign(static_cast<size_t>(NAOB) * NSUB * NPHI, InternalEntry{});

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

KalkerResultCubic KalkerTableCubic::lookup(const double aob, const double ux, const double uy, const double phi) const {
    // Reduce to the stored octant: upsx, upsy >= 0. The mirror relations
    // (reflection about the x-z and y-z planes) map the sign of ux and uy
    // onto the output forces and flip phi (crucially):
    // Fx = sx * Fx_s ,  Fy = sy * Fy_s ,  Mz = sy * Mz_s
    // evaluated at ( |ux|, |uy|, phi_eff = sy * phi ).
    const double sx = (ux < 0) ? -1.0 : 1.0;
    const double sy = (uy < 0) ? -1.0 : 1.0;
    const double axu = std::abs(ux);
    const double ayu = std::abs(uy);

    // ------------------ FIX ------------------
    // Derived from exact symmetries S1 (mirror-y) and S3 (full negation):
    // Fy, Mz are even in ux; Fx is odd in ux.
    // Phi only flips with uy (mirror-y), not ux.
    const double phi_eff = sy * phi;
    // ------------------ FIX ------------------

    const Stencil Sa = make_stencil(pimpl->ax_aob, pimpl->ax_aob.interval(aob));
    const Stencil Sx = make_stencil(pimpl->ax_ups, pimpl->ax_ups.interval(axu));
    const Stencil Sy = make_stencil(pimpl->ax_ups, pimpl->ax_ups.interval(ayu));
    const Stencil Sp = make_stencil(pimpl->ax_phi, pimpl->ax_phi.interval(phi_eff));

    // Successive (tensor-product) monotone-cubic reduction: phi -> uy -> ux -> aob.
    InternalEntry over_a[4];
    for (size_t ia = 0; ia < Sa.n; ++ia) {
        InternalEntry over_x[4];
        for (size_t ix = 0; ix < Sx.n; ++ix) {
            InternalEntry over_y[4];
            for (size_t iy = 0; iy < Sy.n; ++iy) {
                InternalEntry over_p[4];
                for (size_t ip = 0; ip < Sp.n; ++ip) {
                    over_p[ip] = pimpl->at(Sa.col[ia], Sx.col[ix], Sy.col[iy], Sp.col[ip]);
                }
                over_y[iy] = mono_cubic(Sp.x, over_p, Sp.n, phi_eff);
            }
            over_x[ix] = mono_cubic(Sy.x, over_y, Sy.n, ayu);
        }
        over_a[ia] = mono_cubic(Sx.x, over_x, Sx.n, axu);
    }
    const InternalEntry e = mono_cubic(Sa.x, over_a, Sa.n, aob);

    KalkerResultCubic res{};

    // ------------------ FIX ------------------
    res.Fx = sx * e.Fx; // (Fx odd in ux)
    res.Fy = sy * e.Fy; // (Fy even in ux)
    res.Mz = sy * e.Mz; // (Mz even in ux)
    // ------------------ FIX ------------------

    return res;
}
