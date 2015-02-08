/*
 * ============================================================================
 *
 *       Filename:  sim_math.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08.02.2015 00:05:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Bl00drav3n (), 
 *   Organization:  
 *
 * ============================================================================
 */

#include <math.h>

union v2 {

    struct {
        float x;
        float y;
    };

    struct {
        float real;
        float imag;
    };
};

inline v2
V2(float x, float y)
{
    v2 result;
    result.x = x;
    result.y = y;
    return result;
}

inline v2&
operator+=(v2 &a, v2 b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

inline v2&
operator-=(v2 &a, v2 b)
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

inline v2
operator+(v2 a, v2 b)
{
    v2 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    return result;
}

inline v2
operator-(v2 a, v2 b)
{
    v2 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    return result;
}

inline v2
operator*(float f, v2 v)
{
    v2 result;
    result.x = f * v.x;
    result.y = f * v.y;
    return result;
}

inline v2
operator/(v2 v, float f)
{
    v2 result;
    result.x = v.x / f;
    result.y = v.y / f;
    return result;
}

/* utility functions */

inline float
dot(v2 a, v2 b)
{
    float result = a.x * b.x + a.y * b.y;
    return result;
}

inline float
lengthSq(v2 v)
{
    float result = dot(v, v);
    return result;
}

inline float
length(v2 v)
{
    float result;
    result = sqrtf(lengthSq(v));
}

inline float
distanceofSq(v2 a, v2 b)
{
    float result;
    v2 v = a - b;
    result = dot(v, v);
    return result;
}

inline float
distanceof(v2 a, v2 b)
{
    float result;
    result = sqrtf(distanceofSq(a, b));
    return result;
}

/* complex operators */
inline v2
operator*(v2 a, v2 b)
{
    v2 result;
    result.real = a.x * b.x - a.y * b.y;
    result.imag = a.x * b.y + a.y * b.x;
    return result;
}

inline v2&
operator*=(v2 &a, v2 b)
{
    a = a * b;
    return a;
}

inline v2
conjugate(v2 z)
{
    v2 result;
    result = V2(z.real, -z.imag);
    return result;
}

inline v2
operator/(v2 n, v2 d)
{
    v2 result;
    float oneOverDSqr = 1 / lengthSq(d);
    result = oneOverDSqr * n * conjugate(d);
    return result;
}

inline v2
pow(v2 z, uint32_t k)
{
    v2 result = V2(1, 0);
    for(uint32_t i = 0; i < k; i++) {
        result *= z;
    }

    return z;
}

/* general math utilities */

/*
 * TODO: find a better way to calculate those (lookup table?)
 */
static inline uint32_t
binomial(uint32_t n, uint32_t k)
{
    uint32_t result = 1;
    assert(k <= n);
    if(k < n) {
        for(uint32_t i = 1; i <= k; i++) {
            result *= (n + 1 - i) / i;
        }
    }

    return result;
}
