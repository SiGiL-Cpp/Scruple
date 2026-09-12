#include "scruple.h"
#include <cstdio>

using namespace Sigil;

// --- Sanity check: ulp_at reproduces the 750000058 -> 750000064 example ---
// 750000064 sits in [2^29, 2^30), so its ULP is 2^(29-23) = 64.
static_assert(FloatUtils::ulp_at<float>(750000064.f) == 64.f);

// --- next_up / next_down sanity checks -------------------------------------
static_assert(FloatUtils::next_up<float>(0.f) == std::numeric_limits<float>::denorm_min());
static_assert(FloatUtils::next_down<float>(0.f) == -std::numeric_limits<float>::denorm_min());
static_assert(FloatUtils::next_up<float>(1.f) == 1.f + std::numeric_limits<float>::epsilon());
static_assert(FloatUtils::next_down<float>(1.f) == 1.f - std::numeric_limits<float>::epsilon() / 2.f);

// --- The cancellation example from the discussion -----------------------
// Two values pinned at 750000064.f, each declared with twice_abs_err = 64
// (doubled convention -> real bound ±32, i.e. each could really be
// anywhere in [750000032, 750000096]).
using Pinned = Scruple<float, 750000064.f, 750000064.f, 64.f, 0.f>;

constexpr Pinned a{750000064.f};
constexpr Pinned b{750000064.f};
constexpr auto diff = a - b;

// stored difference is exactly 0
static_assert(diff.value() == 0.f);
// Propagated doubled abs_err should be *at least* 128 -> real bound *at
// least* 64, matching the "cumulative absolute error of 64" discussed
// earlier. It's stated as >= rather than == because outward-rounding
// safety margins (needed for soundness in the general case) can add a
// little slack even in a case this clean.
static_assert(decltype(diff)::twice_abs_err >= 128.f);
static_assert(diff.margin() >= 64.f);

// --- Default error at bare declaration -----------------------------------
// This range straddles zero (min_abs == 0) AND reaches into normal
// territory (max_abs == 1500), so it needs BOTH terms simultaneously: a
// nonzero twice_abs_err (the ULP at 0, i.e. denorm_min) to cover the floor
// near zero, and a nonzero twice_rel_err (epsilon) to cover ULP scaling out
// at the large end. One coefficient alone can't cover both regimes.
using Temp = Scruple<float, -1000.f, 1500.f>;
static_assert(Temp::min_abs_value == 0.f);
static_assert(Temp::max_abs_value == 1500.f);
static_assert(Temp::twice_abs_err == std::numeric_limits<float>::denorm_min());
static_assert(Temp::twice_rel_err == std::numeric_limits<float>::epsilon());

// A range that stays entirely in normal, nonzero territory: no floor term
// needed, pure relative error suffices.
using Positive = Scruple<float, 10.f, 20.f>;
static_assert(Positive::min_abs_value == 10.f);
static_assert(Positive::twice_abs_err == 0.f);
static_assert(Positive::twice_rel_err == std::numeric_limits<float>::epsilon());

// A range entirely within subnormal territory: ULP is constant there, so
// only the abs term is needed, rel term correctly suppressed.
using Tiny = Scruple<float, 1e-42f, 1e-41f>;
static_assert(Tiny::twice_rel_err == 0.f);
static_assert(Tiny::twice_abs_err == std::numeric_limits<float>::denorm_min());

// --- Zero-operand soundness fix: an exact-zero error term should not get
// spuriously inflated by outward rounding (this was a real bug: round_up(0)
// used to return denorm_min unconditionally). ---
using ExactZero = Scruple<float, 5.f, 5.f, 0.f, 0.f>;
constexpr ExactZero exact{5.f};
static_assert(exact.margin() == 0.f);

// --- operator+ sanity ------------------------------------------------------
using X = Scruple<float, 0.f, 10.f>;
using Y = Scruple<float, 0.f, 10.f>;
constexpr X x{3.f};
constexpr Y y{4.f};
constexpr auto sum = x + y;
static_assert(decltype(sum)::lower == 0.f);   // 0+0 is exact, no outward slack expected
static_assert(decltype(sum)::upper >= 20.f);  // safe upper bound; outward rounding may add slack
static_assert(sum.value() == 7.f);

// --- default constructor: midpoint, or exactly 0 if the range straddles zero ---
constexpr Scruple<float, 10.f, 20.f> default_mid;
static_assert(default_mid.value() == 15.f);

constexpr Scruple<float, -1000.f, 1500.f> default_straddle;
static_assert(default_straddle.value() == 0.f);

constexpr Scruple<float, 5.f, 5.f, 0.f, 0.f> default_pinned; // Lower==Upper -> midpoint is exactly the pinned value
static_assert(default_pinned.value() == 5.f);

// --- unary minus, and operator- now built entirely on operator+ + unary- ---
constexpr Scruple<float, 2.5f, 1000.f> neg_src{10.f};
constexpr auto negated = -neg_src;
static_assert(decltype(negated)::lower == -1000.f);
static_assert(decltype(negated)::upper == -2.5f);
static_assert(negated.value() == -10.f);

int main() {
    std::printf("diff.value() = %f, margin = +/- %.9g\n", diff.value(), diff.margin());
    std::printf("sum.value()  = %f, range = [%.9g, %.9g]\n",
                sum.value(), decltype(sum)::lower, decltype(sum)::upper);
    std::printf("Temp default: twice_abs_err=%.10g twice_rel_err=%.10g\n",
                Temp::twice_abs_err, Temp::twice_rel_err);
    std::printf("ExactZero margin() = %.10g (should be exactly 0)\n", exact.margin());
    return 0;
}
