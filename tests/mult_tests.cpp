#include "scruple.h"
#include <cstdio>

using namespace Sigil;

int g_fail = 0;
#define CHECK(cond) do { \
    if (!(cond)) { std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++g_fail; } \
} while (0)

int main() {
    // 1. Underflow to zero: both operands' magnitudes so small their exact
    // product (~2.57e-80) is far below the smallest representable
    // subnormal (~1.4e-45), so it rounds to exactly 0.
    {
        using A1 = Scruple<float, 4.6870304754912771e-43f, 4.6870304754912771e-43f>;
        using B1 = Scruple<float, 5.4878151629865804e-38f, 5.4878151629865804e-38f>;
        constexpr A1 a1{4.6870304754912771e-43f};
        constexpr B1 b1{5.4878151629865804e-38f};
        constexpr auto p1 = a1 * b1;
        static_assert(decltype(p1)::lower == 0.f);
        static_assert(decltype(p1)::upper == 0.f);
        CHECK(p1.value() == 0.f);
    }

    // 3. Mixed signs, pinned exact operands: exact product via two_product
    // fast path.
    {
        using A3 = Scruple<float, 1.5f, 1.5f>;
        using B3 = Scruple<float, -2.0f, -2.0f>;
        constexpr A3 a3{1.5f};
        constexpr B3 b3{-2.0f};
        constexpr auto p3 = a3 * b3;
        static_assert(decltype(p3)::lower == -3.f);
        static_assert(decltype(p3)::upper == -3.f);
        CHECK(p3.value() == -3.f);
    }

    // 4. Range crossing zero on one side: 4-corner rule must correctly
    // identify the winning corners (2.0*-3.0=-6 for Min, -1.0*-3.0=3 for
    // Max), not just pair corresponding bounds like +/- does.
    {
        using A4 = Scruple<float, -1.0f, 2.0f>;
        using B4 = Scruple<float, -3.0f, 1.5f>;
        constexpr A4 a4{0.5f};
        constexpr B4 b4{0.0f};
        constexpr auto p4 = a4 * b4;
        static_assert(decltype(p4)::lower == -6.f);
        static_assert(decltype(p4)::upper == 3.f);
    }

    // 7. Denormal * exact one: the "multiply by exact one" fast path must
    // preserve the denormal operand's type EXACTLY unchanged.
    {
        using A7 = Scruple<float, std::numeric_limits<float>::denorm_min(),
                                   std::numeric_limits<float>::denorm_min()>;
        using B7 = Scruple<float, 1.0f, 1.0f>;
        constexpr A7 a7{std::numeric_limits<float>::denorm_min()};
        constexpr B7 b7{1.0f};
        constexpr auto p7 = a7 * b7;
        static_assert(decltype(p7)::lower == std::numeric_limits<float>::denorm_min());
        static_assert(decltype(p7)::upper == std::numeric_limits<float>::denorm_min());
        CHECK(p7.value() == std::numeric_limits<float>::denorm_min());
    }

    // 8. Zero * anything: exact-zero fast path regardless of the other
    // operand's range.
    {
        using A8 = Scruple<float, 0.0f, 0.0f>;
        using B8 = Scruple<float, 42.0f, 42.0f>;
        constexpr A8 a8{0.0f};
        constexpr B8 b8{42.0f};
        constexpr auto p8 = a8 * b8;
        static_assert(decltype(p8)::lower == 0.f);
        static_assert(decltype(p8)::upper == 0.f);
        CHECK(p8.value() == 0.f);
    }

    // 9. Negative underflow: opposite-signed tiny operands, product
    // underflows toward 0 from below.
    {
        using A9 = Scruple<float, -4.6870304754912771e-43f, -4.6870304754912771e-43f>;
        using B9 = Scruple<float, 5.4878151629865804e-38f, 5.4878151629865804e-38f>;
        constexpr A9 a9{-4.6870304754912771e-43f};
        constexpr B9 b9{5.4878151629865804e-38f};
        constexpr auto p9 = a9 * b9;
        static_assert(decltype(p9)::upper == 0.f);
        CHECK(p9.value() == 0.f);
    }

    std::printf("%d failures\n", g_fail);
    return g_fail == 0 ? 0 : 1;
}
