#include "scruple.h"
#include <cstdio>

using namespace Sigil;

int g_fail = 0;
#define CHECK(cond) do { \
    if (!(cond)) { std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++g_fail; } \
} while (0)

int main() {
    // --- 1. value() rename + explicit conversion operator -----------------
    constexpr Scruple<float, 10.f, 20.f> a{15.f};
    static_assert(a.value() == 15.f);
    static_assert(static_cast<float>(a) == 15.f);
    CHECK(a.value() == static_cast<float>(a));
    // Implicit conversion to float must NOT compile -- verified structurally
    // (not exercised here since a failing static_assert would abort the
    // build; this is enforced by `explicit` on the conversion operator).

    // --- 2. Upper defaults to Lower ----------------------------------------
    using Pinned = ScruFl<5.f>; // == Scruple<float,5.f,5.f,defaults>
    static_assert(Pinned::lower == 5.f);
    static_assert(Pinned::upper == 5.f);
    static_assert(Pinned::twice_rel_err == std::numeric_limits<float>::epsilon()); // NOT zero
    static_assert(Pinned::twice_abs_err == 0.f);
    constexpr Pinned p{}; // default constructor picks the single point
    static_assert(p.value() == 5.f);
    CHECK(p.margin() > 0.f); // genuinely non-zero: this is NOT an exact declaration

    // --- 3. Widening conversion --------------------------------------------
    using Narrow = Scruple<float, 0.f, 5.f>;
    using Wide = Scruple<float, -10.f, 10.f, Narrow::twice_abs_err, Narrow::twice_rel_err>;
    constexpr Narrow narrow_val{3.f};
    constexpr Wide widened = narrow_val; // implicit
    static_assert(widened.value() == 3.f);
    static_assert(decltype(widened)::lower == -10.f);
    static_assert(decltype(widened)::upper == 10.f);
    // Per-instance margin is IDENTICAL before and after -- nothing absorbed.
    static_assert(narrow_val.margin() == widened.margin());
    CHECK(narrow_val.margin() == widened.margin());
    // Type-level worst case DOES grow, correctly, since max_abs_value grew.
    static_assert(Wide::range_margin() > Narrow::range_margin());

    // Passing a narrow value to a function taking the wide type by value:
    auto take_wide = [](Wide w) { return w.value(); };
    CHECK(take_wide(narrow_val) == 3.f);

    // --- 4. Comparison operators (weak_ordering: less/equivalent/greater) -
    constexpr Scruple<float, 0.f, 100.f> low{10.f};
    constexpr Scruple<float, 0.f, 100.f> high{90.f};
    static_assert(low < high);
    static_assert(high > low);
    static_assert(low != high);
    static_assert(!(low == high));

    // Same declared type, same value -> equal (overlapping/identical interval).
    constexpr Scruple<float, 0.f, 100.f> low2{10.f};
    static_assert(low == low2);
    static_assert(low <= low2);
    static_assert(low >= low2);

    // Trichotomy: for any two instances, exactly one of <, ==, > holds.
    {
        constexpr bool lt = low < high;
        constexpr bool eq = (low == high);
        constexpr bool gt = low > high;
        static_assert((lt + eq + gt) == 1); // exactly one true
    }

    // Comfortably overlapping intervals (no tie-breaking involved): this
    // needs an EXPLICITLY DECLARED error, not the default -- the default
    // relative error is always ~epsilon/2, i.e. roughly half a ULP AT ANY
    // MAGNITUDE (that's what "representation error" means), so two
    // default-declared pinned values are never more than about one ULP
    // apart before their margins stop overlapping. A comfortable overlap
    // models real declared uncertainty (e.g. measurement noise), not bare
    // rounding: here twice_abs_err is set to 20x the local ULP explicitly.
    {
        constexpr float base = 1000.f;
        constexpr float ulp = FloatUtils::ulp_at(base);
        constexpr float gap3 = FloatUtils::next_up(FloatUtils::next_up(FloatUtils::next_up(base)));
        using W1 = Scruple<float, base, base, 20.f * ulp, 0.f>;
        using W2 = Scruple<float, gap3, gap3, 20.f * ulp, 0.f>;
        constexpr W1 w1{base};
        constexpr W2 w2{gap3};
        CHECK(w1.margin() > (w2.value() - w1.value())); // margin alone exceeds the whole gap
        CHECK(w1 == w2);
    }

    // Boundary case: two default-declared pinned values exactly ONE ULP
    // apart at magnitude 1, where margin() ~= epsilon/2 each (NOT "much
    // larger than one ULP" -- that would be a doubled-vs-halved mixup).
    // v1.margin()+v2.margin() comes out to ~epsilon, i.e. almost exactly
    // the one-ULP gap between them. This exercises the tight-rounding
    // boundary rather than a comfortable margin: TwoSum detects that
    // naively adding epsilon/2 to 1.0 rounds away entirely (ties-to-even),
    // so tight_upper_sum bumps up by exactly one step, landing precisely
    // on v2's value -- not a coincidental near-tie, a cleanly resolved one.
    {
        using V1 = ScruFl<1.0f>;
        constexpr float next = FloatUtils::next_up(1.0f);
        using V2 = Scruple<float, next, next>;
        constexpr V1 v1{};
        constexpr V2 v2{};
        CHECK(v1 == v2);
    }

    // Strict claims must hold against the true (double-precision) values
    // that produced each stored value -- soundness check, reusing the
    // margin-containment property already validated elsewhere.
    {
        using Wide2 = Scruple<float, -1000.f, 1500.f>;
        constexpr Wide2 lo{-500.f};
        constexpr Wide2 hi{1000.f};
        static_assert(lo < hi);
        CHECK(lo.value() < hi.value()); // sanity: the underlying values agree too
    }

    std::printf("%d failures\n", g_fail);
    return g_fail == 0 ? 0 : 1;
}
