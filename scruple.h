#pragma once

#include <bit>
#include <cassert>
#include <compare>
#include <concepts>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <version>

#if defined(__cpp_lib_constexpr_cmath) && __cpp_lib_constexpr_cmath >= 202202L
#include <cmath>
#define SCRUPLE_HAS_CONSTEXPR_FMA 1
#else
#define SCRUPLE_HAS_CONSTEXPR_FMA 0
#endif

namespace Sigil
{
    namespace detail
    {

        template <typename T>
        struct FloatTraits;

        template <>
        struct FloatTraits<float>
        {
            using Bits = std::uint32_t;
            static constexpr int mantissa_bits = 23;
            static constexpr int exponent_bits = 8;
            static constexpr int bias = 127;
        };

        template <>
        struct FloatTraits<double>
        {
            using Bits = std::uint64_t;
            static constexpr int mantissa_bits = 52;
            static constexpr int exponent_bits = 11;
            static constexpr int bias = 1023;
        };

        template <typename T> constexpr T dekker_splitter();
        template <> constexpr float dekker_splitter<float>() { return 4097.0f; }
        template <> constexpr double dekker_splitter<double>() { return 134217729.0; }

        template <typename T>
        struct SplitResult { T hi; T lo; };

        template <typename T>
        constexpr SplitResult<T> split(T a)
        {
            const T c = dekker_splitter<T>() * a;
            const T hi = c - (c - a);
            const T lo = a - hi;
            return {hi, lo};
        }
    } // detail

    namespace FloatUtils
    {
        constexpr auto abs_val(auto x) { return x < decltype(x)(0) ? -x : x; }

        template <typename T>
        constexpr T next_up(T x)
        {
            using Bits = typename detail::FloatTraits<T>::Bits;
            constexpr Bits sign_mask = Bits(1) << (sizeof(Bits) * 8 - 1);
            if (x == T(0)) return std::bit_cast<T>(Bits(1));
            const Bits bits = std::bit_cast<Bits>(x);
            return std::bit_cast<T>((bits & sign_mask) == 0 ? bits + 1 : bits - 1);
        }

        template <typename T>
        constexpr T next_down(T x)
        {
            using Bits = typename detail::FloatTraits<T>::Bits;
            constexpr Bits sign_mask = Bits(1) << (sizeof(Bits) * 8 - 1);
            if (x == T(0)) return std::bit_cast<T>(sign_mask | Bits(1));
            const Bits bits = std::bit_cast<Bits>(x);
            return std::bit_cast<T>((bits & sign_mask) == 0 ? bits - 1 : bits + 1);
        }

        template <typename T>
        constexpr T ulp_at(T x)
        {
            using Traits = detail::FloatTraits<T>;
            using Bits = typename Traits::Bits;

            x = abs_val(x);
            const Bits bits = std::bit_cast<Bits>(x);
            const Bits exp_mask = (Bits(1) << Traits::exponent_bits) - 1;
            const Bits exp_field = (bits >> Traits::mantissa_bits) & exp_mask;

            const Bits smallest_subnormal = Bits(1);

            if (exp_field == 0)
            {
                return std::bit_cast<T>(smallest_subnormal);
            }

            const int unbiased_exp = static_cast<int>(exp_field) - Traits::bias;
            const int ulp_exp = unbiased_exp - Traits::mantissa_bits;
            const int result_exp_field = ulp_exp + Traits::bias;

            if (result_exp_field <= 0)
            {
                return std::bit_cast<T>(smallest_subnormal);
            }

            const Bits result_bits = static_cast<Bits>(result_exp_field) << Traits::mantissa_bits;
            return std::bit_cast<T>(result_bits);
        }

        template <typename T>
        struct TwoSumResult { T sum; T err; };

        template <typename T>
        constexpr TwoSumResult<T> two_sum(T a, T b)
        {
            const T s = a + b;
            const T a_prime = s - b;
            const T b_prime = s - a_prime;
            const T delta_a = a - a_prime;
            const T delta_b = b - b_prime;
            const T err = delta_a + delta_b;
            return {s, err};
        }

        template <typename T>
        struct TwoProductResult { T prod; T err; };

        template <typename T>
        constexpr TwoProductResult<T> two_product(T a, T b)
        {
#if SCRUPLE_HAS_CONSTEXPR_FMA
            const T prod = a * b;
            const T err = std::fma(a, b, -prod);
            return {prod, err};
#else
            const T prod = a * b;
            const auto sa = detail::split(a);
            const auto sb = detail::split(b);
            const T err = ((sa.hi * sb.hi - prod) + sa.hi * sb.lo + sa.lo * sb.hi) + sa.lo * sb.lo;
            return {prod, err};
#endif
        }

        template <typename T> constexpr T tight_upper_sum(T a, T b)
        {
            const auto r = two_sum(a, b);
            return r.err > T(0) ? next_up(r.sum) : r.sum;
        }
        template <typename T> constexpr T tight_lower_sum(T a, T b)
        {
            const auto r = two_sum(a, b);
            return r.err < T(0) ? next_down(r.sum) : r.sum;
        }
        template <typename T> constexpr T tight_upper_diff(T a, T b) { return tight_upper_sum(a, -b); }
        template <typename T> constexpr T tight_lower_diff(T a, T b) { return tight_lower_sum(a, -b); }

        template <typename T> constexpr T tight_upper_product(T a, T b)
        {
            const auto r = two_product(a, b);
            return r.err > T(0) ? next_up(r.prod) : r.prod;
        }
        template <typename T> constexpr T tight_lower_product(T a, T b)
        {
            const auto r = two_product(a, b);
            return r.err < T(0) ? next_down(r.prod) : r.prod;
        }

        template <typename T> constexpr T tight_upper_half(T total)
        {
            if (total == T(0)) return T(0);
            const T half = total / T(2);
            const T reconstructed = half * T(2);
            return reconstructed < total ? next_up(half) : half;
        }
        template <typename T> constexpr T tight_lower_half(T total)
        {
            if (total == T(0)) return T(0);
            const T half = total / T(2);
            const T reconstructed = half * T(2);
            return reconstructed > total ? next_down(half) : half;
        }

        template <typename T> constexpr T tight_upper_quotient(T a, T b)
        {
            const T q = a / b;
            const auto pr = two_product(q, b);
            const T d = pr.prod - a;
            const auto s = two_sum(d, pr.err);
            const T indicator = (s.sum != T(0)) ? s.sum : s.err;
            if (b > T(0))
            {
                return indicator < T(0) ? next_up(q) : q;
            }
            return indicator > T(0) ? next_up(q) : q;
        }

        template <typename T> constexpr T tight_lower_quotient(T a, T b)
        {
            const T q = a / b;
            const auto pr = two_product(q, b);
            const T d = pr.prod - a;
            const auto s = two_sum(d, pr.err);
            const T indicator = (s.sum != T(0)) ? s.sum : s.err;
            if (b > T(0))
            {
                return indicator > T(0) ? next_down(q) : q;
            }
            return indicator < T(0) ? next_down(q) : q;
        }

        template <typename T> constexpr T tight_upper_reciprocal(T x) { return tight_upper_quotient<T>(T(1), x); }
        template <typename T> constexpr T tight_lower_reciprocal(T x) { return tight_lower_quotient<T>(T(1), x); }
    } // FloatUtils

    namespace detail
    {
        template <typename T>
        constexpr T min_abs(T lower, T upper)
        {
            if (lower <= T(0) && upper >= T(0)) return T(0);
            const T al = FloatUtils::abs_val(lower);
            const T au = FloatUtils::abs_val(upper);
            return al < au ? al : au;
        }

        template <typename T>
        constexpr T max_abs(T lower, T upper)
        {
            const T al = FloatUtils::abs_val(lower);
            const T au = FloatUtils::abs_val(upper);
            return al > au ? al : au;
        }

        template <typename T>
        constexpr bool touches_subnormal_or_zero(T lower, T upper)
        {
            return min_abs(lower, upper) < std::numeric_limits<T>::min();
        }

        template <typename T>
        constexpr bool entirely_subnormal_or_zero(T lower, T upper)
        {
            return max_abs(lower, upper) < std::numeric_limits<T>::min();
        }

        template <typename T>
        constexpr T min4(T a, T b, T c, T d)
        {
            T m = a < b ? a : b;
            m = c < m ? c : m;
            m = d < m ? d : m;
            return m;
        }

        template <typename T>
        constexpr T max4(T a, T b, T c, T d)
        {
            T m = a > b ? a : b;
            m = c > m ? c : m;
            m = d > m ? d : m;
            return m;
        }

        template <typename T, T Lower, T Upper>
        constexpr T default_twice_abs_err()
        {
            return touches_subnormal_or_zero<T>(Lower, Upper) ? FloatUtils::ulp_at<T>(min_abs<T>(Lower, Upper)) : T(0);
        }

        template <typename T, T Lower, T Upper>
        constexpr T default_twice_rel_err()
        {
            return entirely_subnormal_or_zero<T>(Lower, Upper) ? T(0) : std::numeric_limits<T>::epsilon();
        }
    }

    template <typename T>
    concept floating_point =
        std::numeric_limits<T>::is_specialized &&
        requires
        {
            typename detail::FloatTraits<T>::Bits;
            { detail::FloatTraits<T>::mantissa_bits } -> std::convertible_to<int>;
            { detail::FloatTraits<T>::exponent_bits } -> std::convertible_to<int>;
            { detail::FloatTraits<T>::bias } -> std::convertible_to<int>;
        } &&
        (sizeof(T) == sizeof(typename detail::FloatTraits<T>::Bits)) &&
        std::regular<T> &&
        std::totally_ordered<T> &&
        requires (T a, T b)
        {
            { a + b } -> std::same_as<T>;
            { a - b } -> std::same_as<T>;
            { a * b } -> std::same_as<T>;
            { a / b } -> std::same_as<T>;
            { -a } -> std::same_as<T>;
        };

    template <typename T, T Lower, T Upper, T TwiceAbsErr, T TwiceRelErr>
    concept ValidScrupleParams =
        floating_point<T> &&
        (Lower <= Upper) &&
        (TwiceAbsErr >= T(0)) &&
        (TwiceRelErr >= T(0));

    template <typename T, T SubLower, T SubUpper, T SupLower, T SupUpper>
    concept RangeSubsetOf =
        (SubLower >= SupLower) &&
        (SubUpper <= SupUpper);

    template <typename T, T Lower, T Upper>
    concept ExcludesZero =
        (Lower > T(0)) || (Upper < T(0));

    template <typename T, T Lower, T Upper = Lower,
            T TwiceAbsErr = detail::default_twice_abs_err<T, Lower, Upper>(),
            T TwiceRelErr = detail::default_twice_rel_err<T, Lower, Upper>()>
        requires ValidScrupleParams<T, Lower, Upper, TwiceAbsErr, TwiceRelErr>
    class Scruple
    {
    public:
        using value_type = T;

        static constexpr T lower = Lower;
        static constexpr T upper = Upper;
        static constexpr T twice_abs_err = TwiceAbsErr;
        static constexpr T twice_rel_err = TwiceRelErr;
        static constexpr T min_abs_value = detail::min_abs(Lower, Upper);
        static constexpr T max_abs_value = detail::max_abs(Lower, Upper);

        constexpr Scruple() noexcept
            : m_value((Lower < T(0) && Upper > T(0)) ? T(0) : (Lower + Upper) / T(2))
        {
        }

        constexpr explicit Scruple(T value) : m_value(value)
        {
            assert(value >= Lower && value <= Upper);
        }

        template <T OtherLower, T OtherUpper>
            requires RangeSubsetOf<T, OtherLower, OtherUpper, Lower, Upper>
        constexpr Scruple(const Scruple<T, OtherLower, OtherUpper, TwiceAbsErr, TwiceRelErr>& other) noexcept
            : m_value(other.value())
        {
        }

        constexpr T value() const { return m_value; }

        constexpr explicit operator T() const noexcept { return m_value; }

        constexpr T margin() const
        {
            const T scaled = FloatUtils::tight_upper_product<T>(twice_rel_err, FloatUtils::abs_val(m_value));
            const T total = FloatUtils::tight_upper_sum<T>(twice_abs_err, scaled);
            return FloatUtils::tight_upper_half<T>(total);
        }

        static constexpr T range_margin()
        {
            const T scaled = FloatUtils::tight_upper_product<T>(TwiceRelErr, max_abs_value);
            const T total = FloatUtils::tight_upper_sum<T>(TwiceAbsErr, scaled);
            return FloatUtils::tight_upper_half<T>(total);
        }

    private:
        T m_value;
    };

    template <typename T, T L1, T U1, T A1, T R1, T L2, T U2, T A2, T R2>
    constexpr auto operator+(const Scruple<T, L1, U1, A1, R1>& lhs,
                            const Scruple<T, L2, U2, A2, R2>& rhs)
    {
        if constexpr (L2 == T(0) && U2 == T(0) && A2 == T(0) && R2 == T(0))
        {
            return Scruple<T, L1, U1, A1, R1>(lhs.value() + rhs.value());
        }
        else if constexpr (L1 == T(0) && U1 == T(0) && A1 == T(0) && R1 == T(0))
        {
            return Scruple<T, L2, U2, A2, R2>(lhs.value() + rhs.value());
        }
        else if constexpr (L1 == U1 && A1 == T(0) && R1 == T(0) &&
                            L2 == U2 && A2 == T(0) && R2 == T(0))
        {
            constexpr auto ts = FloatUtils::two_sum<T>(L1, L2);
            constexpr T new_twice_abs_err = T(2) * FloatUtils::abs_val(ts.err);
            return Scruple<T, ts.sum, ts.sum, new_twice_abs_err, T(0)>(
                lhs.value() + rhs.value());
        }
        else
        {
            constexpr T new_lower = FloatUtils::tight_lower_sum<T>(L1, L2);
            constexpr T new_upper = FloatUtils::tight_upper_sum<T>(U1, U2);

            constexpr T eb1 = Scruple<T, L1, U1, A1, R1>::range_margin();
            constexpr T eb2 = Scruple<T, L2, U2, A2, R2>::range_margin();
            constexpr T propagated = FloatUtils::tight_upper_sum<T>(eb1, eb2);

            constexpr T new_max_abs = detail::max_abs(new_lower, new_upper);
            constexpr T rounding = FloatUtils::ulp_at<T>(new_max_abs);

            constexpr T new_twice_abs_err = FloatUtils::tight_upper_sum<T>(T(2) * propagated, rounding);
            constexpr T new_twice_rel_err = T(0);

            return Scruple<T, new_lower, new_upper, new_twice_abs_err, new_twice_rel_err>(
                lhs.value() + rhs.value());
        }
    }

    template <typename T, T L, T U, T A, T R>
    constexpr Scruple<T, -U, -L, A, R> operator-(const Scruple<T, L, U, A, R>& x)
    {
        return Scruple<T, -U, -L, A, R>(-x.value());
    }

    template <typename T, T L1, T U1, T A1, T R1, T L2, T U2, T A2, T R2>
    constexpr auto operator-(const Scruple<T, L1, U1, A1, R1>& lhs,
                            const Scruple<T, L2, U2, A2, R2>& rhs)
    {
        return lhs + (-rhs);
    }

    template <typename T, T L1, T U1, T A1, T R1, T L2, T U2, T A2, T R2>
    constexpr auto operator*(const Scruple<T, L1, U1, A1, R1>& lhs,
                            const Scruple<T, L2, U2, A2, R2>& rhs)
    {
        if constexpr (L2 == T(0) && U2 == T(0) && A2 == T(0) && R2 == T(0))
        {
            return Scruple<T, T(0), T(0), T(0), T(0)>(lhs.value() * rhs.value());
        }
        else if constexpr (L1 == T(0) && U1 == T(0) && A1 == T(0) && R1 == T(0))
        {
            return Scruple<T, T(0), T(0), T(0), T(0)>(lhs.value() * rhs.value());
        }
        else if constexpr (L2 == T(1) && U2 == T(1) && A2 == T(0) && R2 == T(0))
        {
            return Scruple<T, L1, U1, A1, R1>(lhs.value() * rhs.value());
        }
        else if constexpr (L1 == T(1) && U1 == T(1) && A1 == T(0) && R1 == T(0))
        {
            return Scruple<T, L2, U2, A2, R2>(lhs.value() * rhs.value());
        }
        else if constexpr (L1 == U1 && A1 == T(0) && R1 == T(0) &&
                            L2 == U2 && A2 == T(0) && R2 == T(0))
        {
            constexpr auto tp = FloatUtils::two_product<T>(L1, L2);
            constexpr T new_twice_abs_err = T(2) * FloatUtils::abs_val(tp.err);
            return Scruple<T, tp.prod, tp.prod, new_twice_abs_err, T(0)>(
                lhs.value() * rhs.value());
        }
        else
        {
            constexpr T p_ll_up = FloatUtils::tight_upper_product<T>(L1, L2);
            constexpr T p_lu_up = FloatUtils::tight_upper_product<T>(L1, U2);
            constexpr T p_ul_up = FloatUtils::tight_upper_product<T>(U1, L2);
            constexpr T p_uu_up = FloatUtils::tight_upper_product<T>(U1, U2);
            constexpr T new_upper = detail::max4<T>(p_ll_up, p_lu_up, p_ul_up, p_uu_up);

            constexpr T p_ll_lo = FloatUtils::tight_lower_product<T>(L1, L2);
            constexpr T p_lu_lo = FloatUtils::tight_lower_product<T>(L1, U2);
            constexpr T p_ul_lo = FloatUtils::tight_lower_product<T>(U1, L2);
            constexpr T p_uu_lo = FloatUtils::tight_lower_product<T>(U1, U2);
            constexpr T new_lower = detail::min4<T>(p_ll_lo, p_lu_lo, p_ul_lo, p_uu_lo);

            constexpr T eb1 = Scruple<T, L1, U1, A1, R1>::range_margin();
            constexpr T eb2 = Scruple<T, L2, U2, A2, R2>::range_margin();
            constexpr T max_abs_1 = detail::max_abs(L1, U1);
            constexpr T max_abs_2 = detail::max_abs(L2, U2);

            constexpr T term1 = FloatUtils::tight_upper_product<T>(max_abs_1, eb2);
            constexpr T term2 = FloatUtils::tight_upper_product<T>(max_abs_2, eb1);
            constexpr T term3 = FloatUtils::tight_upper_product<T>(eb1, eb2);
            constexpr T sum12 = FloatUtils::tight_upper_sum<T>(term1, term2);
            constexpr T propagated = FloatUtils::tight_upper_sum<T>(sum12, term3);

            constexpr T new_max_abs = detail::max_abs(new_lower, new_upper);
            constexpr T rounding = FloatUtils::ulp_at<T>(new_max_abs);

            constexpr T new_twice_abs_err = FloatUtils::tight_upper_sum<T>(T(2) * propagated, rounding);
            constexpr T new_twice_rel_err = T(0);

            return Scruple<T, new_lower, new_upper, new_twice_abs_err, new_twice_rel_err>(
                lhs.value() * rhs.value());
        }
    }

    template <typename T, T L, T U, T A, T R>
        requires ExcludesZero<T, L, U>
    constexpr auto inverse(const Scruple<T, L, U, A, R>& x)
    {
        if constexpr (L == U && A == T(0) && R == T(0))
        {
            constexpr T q = T(1) / L;
            constexpr T up_bound = FloatUtils::tight_upper_reciprocal<T>(L);
            constexpr T lo_bound = FloatUtils::tight_lower_reciprocal<T>(L);
            constexpr T new_twice_abs_err = T(2) * (up_bound - lo_bound);
            return Scruple<T, q, q, new_twice_abs_err, T(0)>(T(1) / x.value());
        }
        else
        {
            constexpr T m = detail::min_abs(L, U);
            constexpr T eb_m = FloatUtils::tight_upper_sum<T>(
                A, FloatUtils::tight_upper_product<T>(R, m)) / T(2);
            static_assert(m > eb_m,
                "Scruple::inverse: the range's own declared error reaches too close to "
                "zero for a bounded reciprocal to exist");

            constexpr T new_upper = FloatUtils::tight_upper_reciprocal<T>(L);
            constexpr T new_lower = FloatUtils::tight_lower_reciprocal<T>(U);

            constexpr T denom = FloatUtils::tight_lower_product<T>(m, m - eb_m);
            constexpr T propagated = FloatUtils::tight_upper_quotient<T>(eb_m, denom);

            constexpr T new_max_abs = detail::max_abs(new_lower, new_upper);
            constexpr T rounding = FloatUtils::ulp_at<T>(new_max_abs);

            constexpr T new_twice_abs_err = FloatUtils::tight_upper_sum<T>(T(2) * propagated, rounding);
            constexpr T new_twice_rel_err = T(0);

            return Scruple<T, new_lower, new_upper, new_twice_abs_err, new_twice_rel_err>(
                T(1) / x.value());
        }
    }

    template <typename T, T L1, T U1, T A1, T R1, T L2, T U2, T A2, T R2>
        requires ExcludesZero<T, L2, U2>
    constexpr auto operator/(const Scruple<T, L1, U1, A1, R1>& lhs,
                            const Scruple<T, L2, U2, A2, R2>& rhs)
    {
        return lhs * inverse(rhs);
    }

    template <typename T, T L1, T U1, T A1, T R1, T L2, T U2, T A2, T R2>
    constexpr std::weak_ordering operator<=>(const Scruple<T, L1, U1, A1, R1>& a,
                                            const Scruple<T, L2, U2, A2, R2>& b)
    {
        const T a_hi = FloatUtils::tight_upper_sum<T>(a.value(), a.margin());
        const T b_lo = FloatUtils::tight_lower_diff<T>(b.value(), b.margin());
        if (a_hi < b_lo) return std::weak_ordering::less;

        const T a_lo = FloatUtils::tight_lower_diff<T>(a.value(), a.margin());
        const T b_hi = FloatUtils::tight_upper_sum<T>(b.value(), b.margin());
        if (a_lo > b_hi) return std::weak_ordering::greater;

        return std::weak_ordering::equivalent;
    }

    template <typename T, T L1, T U1, T A1, T R1, T L2, T U2, T A2, T R2>
    constexpr bool operator==(const Scruple<T, L1, U1, A1, R1>& a,
                            const Scruple<T, L2, U2, A2, R2>& b)
    {
        return (a <=> b) == 0;
    }

    template <float Lower, float Upper = Lower,
            float TwiceAbsErr = detail::default_twice_abs_err<float, Lower, Upper>(),
            float TwiceRelErr = detail::default_twice_rel_err<float, Lower, Upper>()>
    using ScruFl = Scruple<float, Lower, Upper, TwiceAbsErr, TwiceRelErr>;

    template <double Lower, double Upper = Lower,
            double TwiceAbsErr = detail::default_twice_abs_err<double, Lower, Upper>(),
            double TwiceRelErr = detail::default_twice_rel_err<double, Lower, Upper>()>
    using ScruDbl = Scruple<double, Lower, Upper, TwiceAbsErr, TwiceRelErr>;

} // Sigil
