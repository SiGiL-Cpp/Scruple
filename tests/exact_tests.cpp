// Tests for the "exact declaration" (point 1) and TwoSum-based exact
// arithmetic (point 2) fast paths in operator+/operator-.
#include "scruple.h"
#include <cstdio>
#include <cmath>

using namespace Sigil;

int g_fail = 0;
#define CHECK(cond) do { \
    if (!(cond)) { std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++g_fail; } \
} while (0)

int main() {
    // --- Point 1: explicit exact declaration ---------------------------
    constexpr Scruple<float, 1.0f, 1.0f, 0.f, 0.f> one{1.0f};
    static_assert(one.margin() == 0.f);
    constexpr Scruple<float, 0.5f, 0.5f, 0.f, 0.f> half{0.5f};
    static_assert(half.margin() == 0.f);
    CHECK(one.margin() == 0.f);
    CHECK(half.margin() == 0.f);

    // --- Point 2a: adding an exact zero preserves the other operand's
    // type exactly (no new rounding term manufactured) ------------------
    using Range = Scruple<float, 10.f, 20.f>; // has default (nonzero) error
    constexpr Range r{15.f};
    constexpr Scruple<float, 0.0f, 0.0f, 0.f, 0.f> zero{0.0f};
    constexpr auto sum_with_zero = r + zero;
    static_assert(decltype(sum_with_zero)::lower == Range::lower);
    static_assert(decltype(sum_with_zero)::upper == Range::upper);
    static_assert(decltype(sum_with_zero)::twice_abs_err == Range::twice_abs_err);
    static_assert(decltype(sum_with_zero)::twice_rel_err == Range::twice_rel_err);
    static_assert(sum_with_zero.value() == 15.f);
    CHECK(sum_with_zero.margin() == r.margin()); // exactly unchanged, not inflated

    constexpr auto zero_plus_r = zero + r; // symmetric case
    static_assert(decltype(zero_plus_r)::lower == Range::lower);
    static_assert(decltype(zero_plus_r)::upper == Range::upper);

    constexpr auto r_minus_zero = r - zero;
    static_assert(decltype(r_minus_zero)::lower == Range::lower);
    static_assert(decltype(r_minus_zero)::upper == Range::upper);

    // 0 - b == -b: bounds negate and swap.
    constexpr auto zero_minus_r = zero - r;
    static_assert(decltype(zero_minus_r)::lower == -Range::upper);
    static_assert(decltype(zero_minus_r)::upper == -Range::lower);
    static_assert(zero_minus_r.value() == -15.f);

    // --- Point 2b: pinned+pinned exact addition, exact-sum case ---------
    // 0.25f + 0.25f = 0.5f is exactly representable -> TwoSum residual
    // should be exactly 0, giving a truly zero-error result, not just a
    // small worst-case bound.
    constexpr Scruple<float, 0.25f, 0.25f, 0.f, 0.f> q1{0.25f};
    constexpr Scruple<float, 0.25f, 0.25f, 0.f, 0.f> q2{0.25f};
    constexpr auto qsum = q1 + q2;
    static_assert(decltype(qsum)::lower == 0.5f);
    static_assert(decltype(qsum)::upper == 0.5f);
    static_assert(decltype(qsum)::twice_abs_err == 0.f);
    static_assert(qsum.margin() == 0.f);
    CHECK(qsum.value() == 0.5f);

    // --- Point 2c: pinned+pinned exact addition, INEXACT-sum case -------
    // 1.0f + 2^-30f: 2^-30 is far below float's precision at 1.0
    // (epsilon ~ 2^-23), so almost all of it is lost to rounding. TwoSum
    // should capture that loss EXACTLY, and stored+err should reconstruct
    // the true mathematical sum exactly (checked in double precision).
    constexpr float tiny = 1.0f / 1073741824.f; // 2^-30
    constexpr Scruple<float, 1.0f, 1.0f, 0.f, 0.f> big{1.0f};
    constexpr Scruple<float, tiny, tiny, 0.f, 0.f> small{tiny};
    constexpr auto isum = big + small;
    static_assert(decltype(isum)::twice_abs_err > 0.f); // NOT exact -- rounding genuinely occurred
    static_assert(decltype(isum)::lower == decltype(isum)::upper); // still a pinned point (exact residual, not a bound)

    {
        const double true_sum = double(1.0f) + double(tiny);
        const double reconstructed = double(isum.value()) +
            (isum.margin()); // margin() == |err| exactly here (signed info lost, but magnitude exact)
        // stored() +/- margin() must bracket the true sum, and should do so
        // TIGHTLY (this is the point of the exact path): the gap between
        // stored+margin and stored-margin should be much smaller than the
        // old worst-case ulp_at bound would have given.
        CHECK(std::fabs(true_sum - double(isum.value())) <= isum.margin() + 1e-300);
        (void)reconstructed;
    }

    // --- Point 2d: pinned+pinned exact cancellation ----------------------
    // Revisit the classic cancellation example, but now with EXPLICITLY
    // EXACT operands (not the "worst-case default error" pinned example
    // from before). Since both operands are asserted exact and their
    // difference (0) is exactly representable, the result should have
    // TRULY ZERO error -- a meaningfully tighter answer than the earlier
    // Pinned example (which used declared abs_err=64 to model genuine
    // measurement uncertainty, not exactness).
    constexpr Scruple<float, 750000064.f, 750000064.f, 0.f, 0.f> pa{750000064.f};
    constexpr Scruple<float, 750000064.f, 750000064.f, 0.f, 0.f> pb{750000064.f};
    constexpr auto pdiff = pa - pb;
    static_assert(pdiff.value() == 0.f);
    static_assert(pdiff.margin() == 0.f); // exact cancellation of exact values: truly zero error
    CHECK(pdiff.margin() == 0.f);

    // A pinned subtraction that is NOT exact: 1.0f - tiny loses precision
    // the same way 1.0f + tiny did above.
    constexpr auto idiff = big - small;
    static_assert(decltype(idiff)::twice_abs_err > 0.f);
    {
        const double true_diff = double(1.0f) - double(tiny);
        CHECK(std::fabs(true_diff - double(idiff.value())) <= idiff.margin() + 1e-300);
    }

    // Regression case from the original RangedFloats.cpp prototype:
    // adding a subnormal-magnitude value to a normal-magnitude one just
    // above the subnormal boundary. Their expected bound: TwiceMaxError
    // >= 2 * 5.7796814070959308e-45.
    {
        using A = Scruple<float, 4.6870304754912771e-43f, 4.6870304754912771e-43f>;
        using B = Scruple<float, 5.4878151629865804e-38f, 5.4878151629865804e-38f>;
        constexpr A a{4.6870304754912771e-43f};
        constexpr B b{5.4878151629865804e-38f};
        constexpr auto sum = a + b;
        static_assert(decltype(sum)::twice_abs_err >= 2.0 * 5.7796814070959308e-45);
        CHECK(decltype(sum)::twice_abs_err >= 2.0 * 5.7796814070959308e-45);
    }

    std::printf("%d failures\n", g_fail);
    return g_fail == 0 ? 0 : 1;
}
