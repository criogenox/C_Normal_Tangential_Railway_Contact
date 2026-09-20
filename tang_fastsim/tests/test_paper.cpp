#include <array>
#include <chrono>
#include <memory>
#include <ostream>
#include <gtest/gtest.h>

#include "fastsim.h"
#include "timing.h"

namespace {
    // Testing data for validation
    struct TestCase {
        Inputs inputs;
        double tol;
        std::array<double, 2> expected;
    };

    // GoogleTest Argument-Dependent Lookup (ADL) matching for output parameter type.
    void PrintTo(const TestCase &tc, std::ostream *os) {
        *os << "{Fx=" << tc.expected[0] << " N, Fy=" << tc.expected[1] << " N}";
    }
}

// 1982 -Kalker- A FAST ALGORITHM FOR THE SIMPLIFIED THEORY OF ROLLING CONTACT (TESTS)
// -----------------------------------------------------------------------------------
// N°0: ux=0 | uy=-2 | fx=2 | fy=4 | mx=5 | my=5 | TOL=3    |=> Tx=0      & Ty=0.1843
// N°1: ux=0 | uy=-2 | fx=2 | fy=4 | mx=5 | my=5 | TOL=0.09 |=> Tx=0      & Ty=0.1850
// N°2: ux=1 | uy=-2 | fx=2 | fy=4 | mx=5 | my=5 | TOL=0.09 |=> Tx=0.5720 & Ty=0.1838
// N°3: ux=1 | uy=-2 | fx=2 | fy=4 | mx=5 | my=5 | TOL=3    |=> Tx=0.5684 & Ty=0.1857
// N°4: ux=1 | uy=-2 | fx=2 | fy=1 | mx=5 | my=5 | TOL=0.09 |=> Tx=0.4530 & Ty=0.8299
// -----------------------------------------------------------------------------------

namespace {
    // ---------------------------------------------------------------------------
    // Parameterized test fixture with pre-warm
    // ---------------------------------------------------------------------------
    class PaperValidation : public ::testing::TestWithParam<TestCase> {
    protected:
        static void SetUpTestSuite() {
            // Warm CPU caches, branch predictor, frequency scaling, allocator.
            // Absorbs one-time costs so the first measured tests isn't penalized.
            const auto p = std::make_unique<Fastsim>(
                Inputs{0.0, -2.0, 2.0, 4.0, 5.0, 5.0, 0.0, 0.0});
            p->creep(3.0);
        }
    };
}

TEST_P(PaperValidation, TestCase) {
    // 0.15 % relative error tolerance on each force component
    constexpr double rel_tol = 0.0015;
    const auto &[inputs, tol, expected] = GetParam();

    const auto p = std::make_unique<Fastsim>(inputs);

    const auto t0 = fastsim::timing::tick_now();
    p->creep(tol);
    const auto t1 = fastsim::timing::tick_now();

    fastsim::timing::bench(t0, t1);

    const auto *info = ::testing::UnitTest::GetInstance()->current_test_info();
    std::cout << "\033[0;38;5;12m* PAPER ("
            << rel_tol * 100.0 << "% tolerance): "
            << info->name() << std::endl;

    // Note Tx sign flip for consistent convention.
    // Definition concern (paper), no algorithm/physics issues.
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const double calculated = (i == 0) ? -p->T[i] : p->T[i];
        ASSERT_NEAR(expected[i], calculated, rel_tol);
    }
}

INSTANTIATE_TEST_SUITE_P(FastSimAlgorithm, PaperValidation,
                         ::testing::Values(
                             TestCase{{0.0, -2.0, 2.0, 4.0, 5.0, 5.0, 0.0, 0.0}, 3.0, {0.0, 0.1843}},
                             TestCase{{0.0, -2.0, 2.0, 4.0, 5.0, 5.0, 0.0, 0.0}, 0.09, {0.0, 0.1850}},
                             TestCase{{1.0, -2.0, 2.0, 4.0, 5.0, 5.0, 0.0, 0.0}, 0.09, {0.5720, 0.1838}},
                             TestCase{{1.0, -2.0, 2.0, 4.0, 5.0, 5.0, 0.0, 0.0}, 3.0, {0.5684, 0.1857}},
                             TestCase{{1.0, -2.0, 2.0, 1.0, 5.0, 5.0, 0.0, 0.0}, 0.09, {0.453, 0.8299}}
                         ));

namespace {
    class FastSimBanner : public ::testing::Environment {
    public:
        void SetUp() override {
            std::cout << "\033[1;38;5;201m# (C++) FASTSIM ALGORITHM FOR ROLLING CONTACT" << std::endl;
            std::cout << "\033[2;38;5;201m# ===========================================" << std::endl;
        }
    };
}

// Registered before main() runs; gtest_main calls RUN_ALL_TESTS().
::testing::Environment *const env =
        ::testing::AddGlobalTestEnvironment(new FastSimBanner);
