/**
 *  @file      dataTypes.h
 *  @brief     Number datatypes Int and Real.
 *  @author    Mikica Kocic
 *  @copyright GNU General Public License (GPLv3).
 */

#ifndef _DATA_TYPES_H_INCLUDED
#define _DATA_TYPES_H_INCLUDED

/////////////////////////////////////////////////////////////////////////////////////////
/** @defgroup g1 Number datatypes                                                      */
/////////////////////////////////////////////////////////////////////////////////////////
                                                                                   /*@{*/
typedef long Int;

/////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_QUADMATH // 128-bit floating point

    #include <quadmath.h>

    typedef __float128  Real;

    #define RealC(v)    v##q
    #define strToReal   strtoflt128
    #define Power(v,e)  powq(v,e)
    #define Sqrt(v)     sqrtq(v)
    #define sqrt(v)     sqrtq(v)
    #define round(v)    roundq(v)

    /** @warning Without abs() is a GCC built-in function, which returns an int.
     */
    inline static Real abs( Real v ) { return fabsq(v); }

    namespace std
    {
        inline static Real isnan( Real v ) { return isnanq(v); }
        inline static Real max( Real a, Real b ) { return fmaxq(a,b); }
    }

    inline static std::ostream& operator<< ( std::ostream& os, const Real& value ) {
        char str[256];
        quadmath_snprintf( str, sizeof(str), "%Qg", value );
        return os << str;
    }

    inline static void fputReal( FILE* outf, Real value ) {
        char str[256];
        quadmath_snprintf( str, 256, "%Qf", value );
        fputs( str, outf );
    }

    inline static Real operator ""_q ( const char* str ) {
        return strtoflt128( str, NULL );
    }

#else
    typedef double      Real;

    #define RealC(v)    v
    #define strToReal   strtod
    #define Power(v,e)  pow(v,e)
    #define Sqrt(v)     sqrt(v)
    #define Erf(v)      erf(v)
    #define Tanh(v)     tanh(v)

    /** @warning Without std::, abs() is a GCC built-in function, which returns an int.
     */
    inline static Real abs( Real v ) { return std::fabs(v); }

    inline static void fputReal( FILE* outf, Real value ) {
        fprintf( outf, "%lf", value );
    }

    inline static Real operator ""_q ( const char* str ) {
        return strtod( str, NULL );
    }
#endif


const Real TINY_Real = 1e-100; // std::numeric_limits<Real>::min();


/** Reads an array of real numbers from a string.
 */
static size_t sscanf_Real( const char* line, Real* data, size_t count )
{
    for( size_t i = 0; i < count; ++ i )
    {
        char* endPtr = NULL;
        data[i] = strToReal( line, &endPtr );
        if( endPtr == line ) return i;
        line = endPtr;
    }

    return count;
}

// static const Real PI = acos(-1.0);

static inline Real pow2( Real v ) { return v * v; }
static inline Real pow3( Real v ) { return v * v * v; }
static inline Real pow4( Real v ) { return v * v * v * v; }

/** Keeps the value of a variable in the given interval.
 */
template <class T> const T& clip( const T& v, const T& minValue, const T& maxValue )
{
    return v < minValue ? minValue : v > maxValue ? maxValue : v;
}

/** Swaps two values.
 */
template <class T> void swap( T& x, T& y ) {
	T temp = x;  x = y;  y = temp;
}
                                                                                   /*@}*/
/////////////////////////////////////////////////////////////////////////////////////////

#endif // _DATA_TYPES_H_INCLUDED
