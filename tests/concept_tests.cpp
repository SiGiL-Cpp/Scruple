#include "scruple.h"
#include <cstdio>

int g_fail = 0;
#define CHECK(cond) do { \
    if (!(cond)) { std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); ++g_fail; } \
} while (0)

// Built-in types Scruple already knows about.
static_assert(Sigil::floating_point<float>);
static_assert(Sigil::floating_point<double>);

// Types with no FloatTraits specialization must be safely rejected --
// SFINAE-style false, not a hard compile error. long double is included
// deliberately: std::floating_point<long double> is true, but our
// concept must still reject it, since we have no FloatTraits<long double>
// specialization (and providing a naive one would be actively wrong, per
// the x86-extended-precision explicit-integer-bit issue discussed earlier).
static_assert(!Sigil::floating_point<int>);
static_assert(!Sigil::floating_point<long double>);
struct NotAFloat {};
static_assert(!Sigil::floating_point<NotAFloat>);

// --- A genuine custom "float-like" type ---------------------------------
// Bit-identical to float (same IEEE754 binary32 layout), just a distinct
// C++ type -- this demonstrates the actual extension mechanism (specialize
// FloatTraits, AND dekker_splitter for the non-FMA fallback path) without
// needing to invent a whole new numeric format.
struct MyFloat
{
    float v;

    constexpr MyFloat() : v(0.0f) {}
    constexpr explicit MyFloat(float x) : v(x) {}

    friend constexpr MyFloat operator+(MyFloat a, MyFloat b) { return MyFloat(a.v + b.v); }
    friend constexpr MyFloat operator-(MyFloat a, MyFloat b) { return MyFloat(a.v - b.v); }
    friend constexpr MyFloat operator*(MyFloat a, MyFloat b) { return MyFloat(a.v * b.v); }
    friend constexpr MyFloat operator/(MyFloat a, MyFloat b) { return MyFloat(a.v / b.v); }
    constexpr MyFloat operator-() const { return MyFloat(-v); }

    friend constexpr bool operator==(MyFloat a, MyFloat b) { return a.v == b.v; }
    friend constexpr auto operator<=>(MyFloat a, MyFloat b) { return a.v <=> b.v; }
};

namespace std
{
    template <>
    struct numeric_limits<MyFloat> : numeric_limits<float>
    {
        static constexpr MyFloat epsilon() noexcept { return MyFloat(numeric_limits<float>::epsilon()); }
        static constexpr MyFloat min() noexcept { return MyFloat(numeric_limits<float>::min()); }
        static constexpr MyFloat denorm_min() noexcept { return MyFloat(numeric_limits<float>::denorm_min()); }
        static constexpr MyFloat lowest() noexcept { return MyFloat(numeric_limits<float>::lowest()); }
        static constexpr MyFloat max() noexcept { return MyFloat(numeric_limits<float>::max()); }
    };
}

namespace Sigil::detail
{
    template <>
    struct FloatTraits<MyFloat>
    {
        using Bits = std::uint32_t;
        static constexpr int mantissa_bits = 23;
        static constexpr int exponent_bits = 8;
        static constexpr int bias = 127;
    };

    template <> constexpr MyFloat dekker_splitter<MyFloat>() { return MyFloat(4097.0f); }
}

// std::bit_cast<Bits>(MyFloat) needs MyFloat to actually BE 32 bits with no
// padding and a trivially-copyable, standard-layout representation matching
// float's bit pattern exactly -- true here since it holds a single float.
static_assert(sizeof(MyFloat) == sizeof(float));

static_assert(Sigil::floating_point<MyFloat>);

int main() {
    using namespace Sigil;

    // Full round trip: declare, add, subtract, compare, using ONLY the
    // custom type -- proving the extension mechanism isn't just a
    // concept-level checkbox, the actual arithmetic works correctly.
    constexpr Scruple<MyFloat, MyFloat(10.0f), MyFloat(20.0f)> a{MyFloat(15.0f)};
    constexpr Scruple<MyFloat, MyFloat(0.0f), MyFloat(10.0f)> b{MyFloat(4.0f)};
    constexpr auto sum = a + b;
    static_assert(sum.value().v == 19.0f);
    static_assert(decltype(sum)::lower.v == 10.0f);
    static_assert(decltype(sum)::upper.v == 30.0f);
    CHECK(sum.value().v == 19.0f);

    constexpr auto neg = -a;
    static_assert(decltype(neg)::lower.v == -20.0f);
    static_assert(decltype(neg)::upper.v == -10.0f);
    CHECK(neg.value().v == -15.0f);

    constexpr Scruple<MyFloat, MyFloat(0.0f), MyFloat(100.0f)> low{MyFloat(10.0f)};
    constexpr Scruple<MyFloat, MyFloat(0.0f), MyFloat(100.0f)> high{MyFloat(90.0f)};
    static_assert(low < high);
    CHECK(low < high);

    // margin() is the thing that actually routes through two_product /
    // dekker_splitter<MyFloat> (via tight_upper_product) -- exercising it
    // here proves the Dekker fallback genuinely runs for a custom type,
    // not just that it compiles unused.
    constexpr auto m = low.margin();
    static_assert(m.v > 0.0f);
    CHECK(m.v > 0.0f);
    std::printf("low.margin() for MyFloat = %.9g (nonzero: two_product path genuinely exercised)\n", m.v);

    std::printf("%d failures\n", g_fail);
    return g_fail == 0 ? 0 : 1;
}
