#pragma once

#include <bit>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace Sigil
{
    namespace detail
    {
        template<typename T>
        strict FindCorrespondingUIntImpl
        {
            using type = std::conditional_t<sizeof(T) == sizeof(uint8_t),  uint8_t,
            std::conditional_t<sizeof(T) == sizeof(uint16_t), uint16_t,
            std::conditional_t<sizeof(T) == sizeof(uint32_t), uint32_t,
            std::conditional_t<sizeof(T) == sizeof(uint64_t), uint64_t,
            #ifdef __SIZEOF_INT128__
            std::conditional_t<sizeof(T) == 16, __uint128_t, void>
            #else
            void
            #endif
            >>>;

            static_assert(!std::is_void<type>, "No matching native unsigned integer type available for this float's size. Customize FloatTraits manually.");
        };

        template<typename T>
        using FindCorrespondingUInt = typename FindCorrespondingUIntImpl<T>::type;

        template <typename T>
        struct FloatTraits
        {
            using Bits = FindCorrespondingUInt<T>;
            static constexpr int mantissa_bits = std::numeric_limits<T>::digits - 1;
            static constexpr int exponent_bits = std::numeric_limits<T>::max_exponent - std::numeric_limits<T>::min_exponent + 2;
            static constexpr int bias = (1 << (exponent_bits - 1)) - 1;
        };

        static_assert(std::is_same_v<FloatTraits<float>::Bits, std::uint32_t>);
        static_assert(FloatTraits<float>::mantissa_bits == 23);
        static_assert(FloatTraits<float>::exponent_bits == 8);
        static_assert(FloatTraits<float>::bias == 127);

        static_assert(std::is_same_v<FloatTraits<double>::Bits, std::uint64_t>);
        static_assert(FloatTraits<double>::mantissa_bits == 52);
        static_assert(FloatTraits<double>::exponent_bits == 11);
        static_assert(FloatTraits<double>::bias == 1023);

        constexpr auto abs_val(auto x) { return x < decltype(x)(0) ? -x : x; }

        template <typename T>
        constexpr T next_up(T x)
        {
            using Bits = typename FloatTraits<T>::Bits;
            constexpr Bits sign_mask = Bits(1) << (sizeof(Bits) * 8 - 1);
            if (x == T(0)) return std::bit_cast<T>(Bits(1));
            const Bits bits = std::bit_cast<Bits>(x);
            return std::bit_cast<T>((bits & sign_mask) == 0 ? bits + 1 : bits - 1);
        }

        template <typename T>
        constexpr T next_down(T x)
        {
            using Bits = typename FloatTraits<T>::Bits;
            constexpr Bits sign_mask = Bits(1) << (sizeof(Bits) * 8 - 1);
            if (x == T(0)) return std::bit_cast<T>(sign_mask | Bits(1));
            const Bits bits = std::bit_cast<Bits>(x);
            return std::bit_cast<T>((bits & sign_mask) == 0 ? bits - 1 : bits + 1);
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

        template <typename T>
        struct TwoProductResult { T prod; T err; };

        template <typename T>
        constexpr TwoProductResult<T> two_product(T a, T b)
        {
            const T prod = a * b;
            const auto sa = split(a);
            const auto sb = split(b);
            const T err = ((sa.hi * sb.hi - prod) + sa.hi * sb.lo + sa.lo * sb.hi) + sa.lo * sb.lo;
            return {prod, err};
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

        template <typename T>
        constexpr T ulp_at(T x)
        {
            using Traits = FloatTraits<T>;
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
        constexpr T min_abs(T lower, T upper)
        {
            if (lower <= T(0) && upper >= T(0)) return T(0);
            const T al = abs_val(lower);
            const T au = abs_val(upper);
            return al < au ? al : au;
        }

        template <typename T>
        constexpr T max_abs(T lower, T upper)
        {
            const T al = abs_val(lower);
            const T au = abs_val(upper);
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

        template <typename T, T Lower, T Upper>
        constexpr T default_twice_abs_err()
        {
            return touches_subnormal_or_zero<T>(Lower, Upper) ? ulp_at<T>(min_abs<T>(Lower, Upper)) : T(0);
        }

        template <typename T, T Lower, T Upper>
        constexpr T default_twice_rel_err()
        {
            return entirely_subnormal_or_zero<T>(Lower, Upper) ? T(0) : std::numeric_limits<T>::epsilon();
        }
    } // detail

    template <typename T, T Lower, T Upper,
            T TwiceAbsErr = detail::default_twice_abs_err<T, Lower, Upper>(),
            T TwiceRelErr = detail::default_twice_rel_err<T, Lower, Upper>()>
    class Scruple
    {
        static_assert(std::is_floating_point_v<T>, "Scruple requires a floating-point type");
        static_assert(Lower <= Upper, "Scruple: lower bound must not exceed upper bound");
        static_assert(TwiceAbsErr >= T(0) && TwiceRelErr >= T(0), "Scruple: error terms must be non-negative");

    public:
        using value_type = T;

        static constexpr T lower = Lower;
        static constexpr T upper = Upper;
        static constexpr T twice_abs_err = TwiceAbsErr;
        static constexpr T twice_rel_err = TwiceRelErr;
        static constexpr T min_abs_value = detail::min_abs(Lower, Upper);
        static constexpr T max_abs_value = detail::max_abs(Lower, Upper);

        constexpr explicit Scruple(T value) : m_value(value) {}

        constexpr T stored() const { return m_value; }

        constexpr T margin() const
        {
            const T scaled = detail::tight_upper_product<T>(twice_rel_err, detail::abs_val(m_value));
            const T total = detail::tight_upper_sum<T>(twice_abs_err, scaled);
            return detail::tight_upper_half<T>(total);
        }

        static constexpr T range_margin()
        {
            const T scaled = detail::tight_upper_product<T>(TwiceRelErr, max_abs_value);
            const T total = detail::tight_upper_sum<T>(TwiceAbsErr, scaled);
            return detail::tight_upper_half<T>(total);
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
            return Scruple<T, L1, U1, A1, R1>(lhs.stored() + rhs.stored());
        }
        else if constexpr (L1 == T(0) && U1 == T(0) && A1 == T(0) && R1 == T(0))
        {
            return Scruple<T, L2, U2, A2, R2>(lhs.stored() + rhs.stored());
        }
        else if constexpr (L1 == U1 && A1 == T(0) && R1 == T(0) &&
                            L2 == U2 && A2 == T(0) && R2 == T(0))
        {
            constexpr auto ts = detail::two_sum<T>(L1, L2);
            constexpr T new_twice_abs_err = T(2) * detail::abs_val(ts.err);
            return Scruple<T, ts.sum, ts.sum, new_twice_abs_err, T(0)>(
                lhs.stored() + rhs.stored());
        }
        else
        {
            constexpr T new_lower = detail::tight_lower_sum<T>(L1, L2);
            constexpr T new_upper = detail::tight_upper_sum<T>(U1, U2);

            constexpr T eb1 = Scruple<T, L1, U1, A1, R1>::range_margin();
            constexpr T eb2 = Scruple<T, L2, U2, A2, R2>::range_margin();
            constexpr T propagated = detail::tight_upper_sum<T>(eb1, eb2);

            constexpr T new_max_abs = detail::max_abs(new_lower, new_upper);
            constexpr T rounding = detail::ulp_at<T>(new_max_abs);

            constexpr T new_twice_abs_err = detail::tight_upper_sum<T>(T(2) * propagated, rounding);
            constexpr T new_twice_rel_err = T(0);

            return Scruple<T, new_lower, new_upper, new_twice_abs_err, new_twice_rel_err>(
                lhs.stored() + rhs.stored());
        }
    }

    template <typename T, T L1, T U1, T A1, T R1, T L2, T U2, T A2, T R2>
    constexpr auto operator-(const Scruple<T, L1, U1, A1, R1>& lhs,
                            const Scruple<T, L2, U2, A2, R2>& rhs)
    {
        if constexpr (L2 == T(0) && U2 == T(0) && A2 == T(0) && R2 == T(0))
        {
            return Scruple<T, L1, U1, A1, R1>(lhs.stored() - rhs.stored());
        }
        else if constexpr (L1 == T(0) && U1 == T(0) && A1 == T(0) && R1 == T(0))
        {
            return Scruple<T, -U2, -L2, A2, R2>(lhs.stored() - rhs.stored());
        }
        else if constexpr (L1 == U1 && A1 == T(0) && R1 == T(0) &&
                            L2 == U2 && A2 == T(0) && R2 == T(0))
        {
            constexpr auto ts = detail::two_sum<T>(L1, -L2);
            constexpr T new_twice_abs_err = T(2) * detail::abs_val(ts.err);
            return Scruple<T, ts.sum, ts.sum, new_twice_abs_err, T(0)>(
                lhs.stored() - rhs.stored());
        }
        else
        {
            constexpr T new_lower = detail::tight_lower_diff<T>(L1, U2);
            constexpr T new_upper = detail::tight_upper_diff<T>(U1, L2);

            constexpr T eb1 = Scruple<T, L1, U1, A1, R1>::range_margin();
            constexpr T eb2 = Scruple<T, L2, U2, A2, R2>::range_margin();
            constexpr T propagated = detail::tight_upper_sum<T>(eb1, eb2);

            constexpr T new_max_abs = detail::max_abs(new_lower, new_upper);
            constexpr T rounding = detail::ulp_at<T>(new_max_abs);

            constexpr T new_twice_abs_err = detail::tight_upper_sum<T>(T(2) * propagated, rounding);
            constexpr T new_twice_rel_err = T(0);

            return Scruple<T, new_lower, new_upper, new_twice_abs_err, new_twice_rel_err>(
                lhs.stored() - rhs.stored());
        }
    }

    template <float Lower, float Upper,
            float TwiceAbsErr = detail::default_twice_abs_err<float, Lower, Upper>(),
            float TwiceRelErr = detail::default_twice_rel_err<float, Lower, Upper>()>
    using ScruFl = Scruple<float, Lower, Upper, TwiceAbsErr, TwiceRelErr>;

    template <double Lower, double Upper,
            double TwiceAbsErr = detail::default_twice_abs_err<double, Lower, Upper>(),
            double TwiceRelErr = detail::default_twice_rel_err<double, Lower, Upper>()>
    using ScruDbl = Scruple<double, Lower, Upper, TwiceAbsErr, TwiceRelErr>;

} // Sigil
