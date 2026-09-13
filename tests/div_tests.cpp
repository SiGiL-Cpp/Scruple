#include "scruple.h"
#include <cstdio>
#include <random>
#include <cmath>

using namespace Sigil;

int g_fail = 0;
#define CHECK(cond) do { \
    if (!(cond)) { std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++g_fail; } \
} while (0)

int main() {
    // --- Pinned-exact inverse: exact via tight_upper/lower_reciprocal ----
    {
        using A = Scruple<float, 4.0f, 4.0f, 0.f, 0.f>;
        constexpr A a{4.0f};
        constexpr auto inv = inverse(a);
        static_assert(decltype(inv)::lower == 0.25f);
        static_assert(decltype(inv)::upper == 0.25f);
        static_assert(inv.margin() == 0.f); // 1/4 is exactly representable
        CHECK(inv.value() == 0.25f);
    }
    {
        // 1/3 is NOT exactly representable -- margin should be nonzero but tiny.
        using A = Scruple<float, 3.0f, 3.0f, 0.f, 0.f>;
        constexpr A a{3.0f};
        constexpr auto inv = inverse(a);
        static_assert(decltype(inv)::twice_abs_err > 0.f);
        CHECK(inv.margin() > 0.f);
        CHECK(inv.margin() < 1e-6f); // should be on the order of 1 ULP at 1/3, not larger
        double true_recip = 1.0 / 3.0;
        CHECK(std::fabs(true_recip - static_cast<double>(inv.value())) <= static_cast<double>(inv.margin()));
    }

    // --- Negative pinned value ---------------------------------------------
    {
        using A = Scruple<float, -2.0f, -2.0f, 0.f, 0.f>;
        constexpr A a{-2.0f};
        constexpr auto inv = inverse(a);
        static_assert(decltype(inv)::lower == -0.5f);
        static_assert(decltype(inv)::upper == -0.5f);
        CHECK(inv.value() == -0.5f);
    }

    // --- Range inverse: order flips correctly ([1,2] -> [0.5,1]) -----------
    {
        using A = Scruple<float, 1.0f, 2.0f, 0.f, 0.f>;
        constexpr A a{1.5f};
        constexpr auto inv = inverse(a);
        static_assert(decltype(inv)::lower == 0.5f);
        static_assert(decltype(inv)::upper == 1.0f);
    }
    // Negative range: [-2,-1] -> [-1,-0.5]
    {
        using A = Scruple<float, -2.0f, -1.0f, 0.f, 0.f>;
        constexpr A a{-1.5f};
        constexpr auto inv = inverse(a);
        static_assert(decltype(inv)::lower == -1.0f);
        static_assert(decltype(inv)::upper == -0.5f);
    }

    // --- operator/ built on inverse() ----------------------------------
    {
        using A = Scruple<float, 6.0f, 6.0f, 0.f, 0.f>;
        using B = Scruple<float, 2.0f, 2.0f, 0.f, 0.f>;
        constexpr A a{6.0f};
        constexpr B b{2.0f};
        constexpr auto q = a / b;
        static_assert(decltype(q)::lower == 3.f);
        static_assert(decltype(q)::upper == 3.f);
        CHECK(q.value() == 3.f);
    }

    // --- Monte Carlo soundness: division against a double-precision
    // ground truth, for a range strictly excluding zero. ------------------
    {
        std::mt19937_64 gen(0xD1D1DE);
        std::uniform_real_distribution<double> da(-100.0, 100.0);
        std::uniform_real_distribution<double> db(10.0, 50.0); // strictly positive divisor range
        int fails = 0;
        int trials = 20000;
        double min_slack = 1.0;

        using B = Scruple<float, 10.0f, 50.0f>;

        for (int i = 0; i < trials; ++i) {
            double ta = da(gen);
            double tb = db(gen);
            float sa = static_cast<float>(ta);
            float sb = static_cast<float>(tb);

            Scruple<float, -100.0f, 100.0f> a(sa);
            B b(sb);
            auto q = a / b;

            double true_q = ta / tb;
            double bound = static_cast<double>(q.margin());
            double actual_err = std::fabs(true_q - static_cast<double>(q.value()));
            double lo = static_cast<double>(decltype(q)::lower);
            double hi = static_cast<double>(decltype(q)::upper);

            if (!(actual_err <= bound)) {
                if (fails < 5) std::printf("UNSOUND div: ta=%.17g tb=%.17g true_q=%.17g stored=%.9g err=%.3e bound=%.3e\n",
                    ta, tb, true_q, q.value(), actual_err, bound);
                ++fails;
            }
            if (!(true_q >= lo && true_q <= hi)) {
                if (fails < 5) std::printf("UNSOUND range: true_q=%.17g not in [%.9g,%.9g]\n", true_q, lo, hi);
                ++fails;
            }
            if (bound > 0.0) min_slack = std::min(min_slack, 1.0 - actual_err / bound);
        }
        std::printf("division soundness: %d fails out of %d trials, min slack ratio=%.4f\n",
                     fails, trials, min_slack);
        CHECK(fails == 0);
    }

    std::printf("%d failures\n", g_fail);
    return g_fail == 0 ? 0 : 1;
}
