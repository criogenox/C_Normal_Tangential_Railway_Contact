#include <array>
#include <chrono>
#include <cmath>
#include <memory>
#include <numbers>
#include <ostream>
#include <gtest/gtest.h>

#include "fastsim.h"
#include "timing.h"

// Physical constants shared across all CONTACT validation cases
// G [Pa-N/m], mu = friction coefficient, nu = Poisson's ratio
// (note: C11/C22/C23 are pre-fitted for nu = 0.27–0.28)
namespace gen {
    inline constexpr double G = 82e9;
    inline constexpr double mu = 0.3;
    inline constexpr double nu = 0.28;
    inline constexpr double a = 0.007257; // [m]
    inline constexpr double b = 0.005155; // [m]
}

namespace {
    // Test case: kinematic creepage parameters + expected forces from CONTACT software
    struct TestCase {
        double gx; // Longitudinal creepage [-]
        double gy; // Lateral creepage [-]
        double wz_permm; // Spin [1/mm]
        double N; // Normal load [N]
        std::array<double, 2> expected; // {Fx, Fy} in [N]
    };

    // GoogleTest Argument-Dependent Lookup (ADL) matching for output parameter type.
    void PrintTo(const TestCase &tc, std::ostream *os) {
        *os << "{Fx=" << tc.expected[0] << " N, Fy=" << tc.expected[1] << " N}";
    }
}

// CONTACT VALIDATION (TESTS) - PHYSICAL PARAMETERS: [-] [-] [rad/mm] [N] [N] [N]
// ------------------------------------------------------------------------------
// N°1: gx=0.002 | gy=0.003 | wz=0.0001 | N=90e3 |=> Fx=-14350.0 & Fy=-22320.0
// N°2: gx=0.03  | gy=0.02  | wz=0.001  | N=90e3 |=> Fx=-22270.0 & Fy=-15110.0
// N°3: gx=0.003 | gy=0.002 | wz=0.0001 | N=45e3 |=> Fx=-10960.0 & Fy=-7827.0
// N°4: gx=0.003 | gy=0.002 | wz=0.0001 | N=90e3 |=> Fx=-22058.0 & Fy=-15190.0
// ------------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Normalize physical parameters → dimensionless FASTSIM inputs
// ---------------------------------------------------------------------------
static Inputs normalize(const TestCase &tc) {
    const double muN = gen::mu * tc.N;
    const double wz = tc.wz_permm * 1000.0; // rad/mm → rad/m

    constexpr double rel = gen::b / gen::a;
    constexpr double ab = gen::a * gen::b;

    // Kalker's creep coefficients (valid for 0.04 < b/a < 25, nu ≈ 0.27–0.28)
    const double C11 = 3.2893 + 0.975 / rel - 0.012 / std::pow(rel, 2);
    const double C22 = 2.4014 + 1.3179 / rel - 0.02 / std::pow(rel, 2);
    const double C23 = 0.4147 + 1.0184 / rel + 0.0565 / std::pow(rel, 2)
                       - 0.0013 / std::pow(rel, 3);

    // Dimensionless creepages
    const double ux = 3.0 * std::numbers::pi * ab * gen::G * C11 * tc.gx / (16.0 * muN);
    const double uy = 3.0 * std::numbers::pi * ab * gen::G * C22 * tc.gy / (16.0 * muN);
    const double fy = 2.0 * std::pow(ab, 1.5) * gen::G * C23 * wz / muN;
    const double fx = fy * (gen::b / gen::a);

    return Inputs{
        .ux = ux, .uy = uy, .fx = fx, .fy = fy,
        .mx = 5.0, .my = 5.0, .Tx = 0.0, .Ty = 0.0
    };
}

namespace {
    // ---------------------------------------------------------------------------
    // Parameterized test fixture with pre-warm
    // ---------------------------------------------------------------------------
    class ContactValidation : public ::testing::TestWithParam<TestCase> {
    protected:
        static void SetUpTestSuite() {
            // Warm CPU caches, branch predictor, frequency scaling, allocator.
            // Absorbs one-time costs so the first measured tests isn't penalized.
            const auto p = std::make_unique<Fastsim>(
                normalize(TestCase{0.002, 0.003, 0.0001, 90e3, {}}));
            p->creep(1e-3);
        }
    };
}

TEST_P(ContactValidation, TestCase) {
    const auto &tc = GetParam();

    const auto inputs = normalize(tc);
    const auto p = std::make_unique<Fastsim>(inputs);

    const auto t0 = fastsim::timing::tick_now();
    p->creep(1e-3);
    const auto t1 = fastsim::timing::tick_now();

    fastsim::timing::bench(t0, t1);

    // Normalized results back to physical forces [N]
    const double muN = gen::mu * tc.N;
    const double Fx = p->T[0] * muN;
    const double Fy = p->T[1] * muN;

    // 3.5 % relative error tolerance on each force component
    constexpr double rel_tol = 0.035;
    const double fx_tol = std::abs(tc.expected[0]) * rel_tol;
    const double fy_tol = std::abs(tc.expected[1]) * rel_tol;

    const auto *info = ::testing::UnitTest::GetInstance()->current_test_info();
    std::cout << "\033[0;38;5;12m* CONTACT ("
            << rel_tol * 100.0 << "% tolerance): "
            << info->name() << std::endl;

    EXPECT_NEAR(tc.expected[0], Fx, fx_tol)
        << "Fx outside 3.5% tolerance";
    EXPECT_NEAR(tc.expected[1], Fy, fy_tol)
        << "Fy outside 3.5% tolerance";
}

INSTANTIATE_TEST_SUITE_P(FastSimAlgorithm, ContactValidation,
                         ::testing::Values(
                             // Input: gx [-] | gy [-] | wz [rad/mm] | N (load [N]) | Expected: Fx [N] | Fy [N]
                             TestCase{0.002, 0.003, 0.0001, 90e3, {-14350.0, -22320.0}},
                             TestCase{0.03, 0.02, 0.001, 90e3, {-22270.0, -15110.0}},
                             TestCase{0.003, 0.002, 0.0001, 45e3, {-10960.0, -7827.0}},
                             TestCase{0.003, 0.002, 0.0001, 90e3, {-22058.0, -15190.0}}
                         ));
