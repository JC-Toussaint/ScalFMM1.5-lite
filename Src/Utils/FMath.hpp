// See LICENCE file at project root
#ifndef FMATH_HPP
#define FMATH_HPP

#include <cmath>
#include <limits>

#include "FGlobal.hpp"

#ifdef SCALFMM_USE_SSE
#include "FSse.hpp"
#endif

#ifdef SCALFMM_USE_AVX
#include "FAvx.hpp"
#endif

/**
 * @author Berenger Bramas (berenger.bramas@inria.fr)
 * @class
 * Please read the license
 *
 * Propose basic math functions or indirections to std math.
 */
struct FMath{
    template <class FReal>
    constexpr static FReal FPi(){ return FReal(M_PI); }
    template <class FReal>
    constexpr static FReal FTwoPi(){ return FReal(2.0*M_PI); }
    template <class FReal>
    constexpr static FReal FPiDiv2(){ return FReal(M_PI_2); }

    /** To get absolute value */
    template <class NumType>
    static NumType Abs(const NumType inV){
        return (inV < 0 ? -inV : inV);
    }

    /** To know if 2 values seems to be equal */

    /** To get pow */
    static double pow(double x, double y){
        return ::pow(x,y);
    }
    static float pow(float x, float y){
        return ::powf(x,y);
    }
    template <class NumType>
    static NumType pow(const NumType inValue, int power){
        NumType result = 1;
        while(power-- > 0) result *= inValue;
        return result;
    }

    /** To know if a value is between two others */
    template <class NumType>
    static bool Between(const NumType inValue, const NumType inMin, const NumType inMax){
        return ( inMin <= inValue && inValue < inMax );
    }
    /** To compute fmadd operations **/
    template <class NumType>
    static NumType FMAdd(const NumType a, const NumType b, const NumType c){
	return a * b + c;
    }

#if  defined(SCALFMM_USE_SSE ) && defined(__SSSE4_1__)
    static __m128 FMAdd(const __m128 inV1, const __m128 inV2, const __m128 inV3){
        return _mm_add_ps( _mm_mul_ps(inV1,inV2), inV3);
    }

    static __m128d FMAdd(const __m128d inV1, const __m128d inV2, const __m128d inV3){
        return _mm_add_pd( _mm_mul_pd(inV1,inV2), inV3);
    }

#endif
#ifdef SCALFMM_USE_AVX
    static __m256 FMAdd(const __m256 inV1, const __m256 inV2, const __m256 inV3){
        return _mm256_add_ps( _mm256_mul_ps(inV1,inV2), inV3);
    }

    static __m256d FMAdd(const __m256d inV1, const __m256d inV2, const __m256d inV3){
        return _mm256_add_pd( _mm256_mul_pd(inV1,inV2), inV3);
    }

#endif
    /** To get sqrt of a FReal */
    static float Sqrt(const float inValue){
        return sqrtf(inValue);
    }
    static double Sqrt(const double inValue){
        return sqrt(inValue);
    }
    static float Rsqrt(const float inValue){
        return float(1.0)/sqrtf(inValue);
    }
    static double Rsqrt(const double inValue){
        return 1.0/sqrt(inValue);
    }
#ifdef SCALFMM_USE_SSE

    static __m128 Sqrt(const __m128 inV){
        return _mm_sqrt_ps(inV);
    }

    static __m128d Sqrt(const __m128d inV){
        return _mm_sqrt_pd(inV);
    }

    static __m128 Rsqrt(const __m128 inV){
        return _mm_rsqrt_ps(inV);
    }

    static __m128d Rsqrt(const __m128d inV){
        return _mm_set_pd1(1.0) / _mm_sqrt_pd(inV);
    }
#endif
#ifdef SCALFMM_USE_AVX

    static __m256 Sqrt(const __m256 inV){
        return _mm256_sqrt_ps(inV);
    }

    static __m256d Sqrt(const __m256d inV){
        return _mm256_sqrt_pd(inV);
    }

    static __m256 Rsqrt(const __m256 inV){
        return _mm256_rsqrt_ps(inV);
    }

    static __m256d Rsqrt(const __m256d inV){
        return _mm256_set1_pd(1.0) / _mm256_sqrt_pd(inV);
    }
#endif

    /** To get atan2 of a 2 FReal,  The return value is given in radians and is in the
      range -pi to pi, inclusive.  */
    static float Atan2(const float inValue1,const float inValue2){
        return atan2f(inValue1,inValue2);
    }
    static double Atan2(const double inValue1,const double inValue2){
        return atan2(inValue1,inValue2);
    }

    template <class NumType>
    static NumType Zero();

    template <class NumType>
    static NumType One();

    template <class DestType, class SrcType>
    static DestType ConvertTo(const SrcType val);

};

template <>
inline float FMath::Zero<float>(){
    return float(0.0);
}

template <>
inline double FMath::Zero<double>(){
    return double(0.0);
}

template <>
inline float FMath::One<float>(){
    return float(1.0);
}

template <>
inline double FMath::One<double>(){
    return double(1.0);
}

template <>
inline float FMath::ConvertTo<float,float>(const float val){
    return val;
}

template <>
inline double FMath::ConvertTo<double,double>(const double val){
    return val;
}

template <>
inline float FMath::ConvertTo<float,const float*>(const float* val){
    return *val;
}

template <>
inline double FMath::ConvertTo<double,const double*>(const double* val){
    return *val;
}

#ifdef SCALFMM_USE_SSE
template <>
inline __m128 FMath::One<__m128>(){
    return _mm_set_ps1(1.0);
}

template <>
inline __m128d FMath::One<__m128d>(){
    return _mm_set_pd1(1.0);
}

template <>
inline __m128 FMath::Zero<__m128>(){
    return _mm_setzero_ps();
}

template <>
inline __m128d FMath::Zero<__m128d>(){
    return _mm_setzero_pd();
}

template <>
inline __m128 FMath::ConvertTo<__m128,float>(const float val){
    return _mm_set_ps1(val);
}

template <>
inline __m128d FMath::ConvertTo<__m128d,double>(const double val){
    return _mm_set_pd1(val);
}

template <>
inline __m128 FMath::ConvertTo<__m128,const float*>(const float* val){
    return _mm_load1_ps(val);
}

template <>
inline __m128d FMath::ConvertTo<__m128d,const double*>(const double* val){
    return _mm_load1_pd(val);
}

template <>
inline float FMath::ConvertTo<float,__m128>(const __m128 val){
    __attribute__((aligned(16))) float buffer[4];
    _mm_store_ps(buffer, val);
    return buffer[0] + buffer[1] + buffer[2] + buffer[3];
}

template <>
inline double FMath::ConvertTo<double,__m128d>(const __m128d val){
    __attribute__((aligned(16))) double buffer[2];
    _mm_store_pd(buffer, val);
    return buffer[0] + buffer[1];
}
#endif

#ifdef SCALFMM_USE_AVX
template <>
inline __m256 FMath::One<__m256>(){
    return _mm256_set1_ps(1.0);
}

template <>
inline __m256d FMath::One<__m256d>(){
    return _mm256_set1_pd(1.0);
}

template <>
inline __m256 FMath::Zero<__m256>(){
    return _mm256_setzero_ps();
}

template <>
inline __m256d FMath::Zero<__m256d>(){
    return _mm256_setzero_pd();
}

template <>
inline __m256 FMath::ConvertTo<__m256,float>(const float val){
    return _mm256_set1_ps(val);
}

template <>
inline __m256d FMath::ConvertTo<__m256d,double>(const double val){
    return _mm256_set1_pd(val);
}

template <>
inline __m256 FMath::ConvertTo<__m256,const float*>(const float* val){
    return _mm256_broadcast_ss(val);
}

template <>
inline __m256d FMath::ConvertTo<__m256d,const double*>(const double* val){
    return _mm256_broadcast_sd(val);
}

template <>
inline float FMath::ConvertTo<float,__m256>(const __m256 val){
    __attribute__((aligned(32))) float buffer[8];
    _mm256_store_ps(buffer, val);
    return buffer[0] + buffer[1] + buffer[2] + buffer[3] + buffer[4] + buffer[5] + buffer[6] + buffer[7];
}

template <>
inline double FMath::ConvertTo<double,__m256d>(const __m256d val){
    __attribute__((aligned(32))) double buffer[4];
    _mm256_store_pd(buffer, val);
    return buffer[0] + buffer[1] + buffer[2] + buffer[3];
}
#endif

#endif //FMATH_HPP
