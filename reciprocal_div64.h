#ifndef RECIPROCAL_DIV64_H
#define RECIPROCAL_DIV64_H

#include <immintrin.h>

#include "bitops.h"

#ifdef WIN32
        #include "int128_c.h"
#endif


/*
 * This algorithm is based on the paper "Division by Invariant
 * Integers Using Multiplication" by Torbjörn Granlund and Peter
 * L. Montgomery.
 *
 * The assembler implementation from Agner Fog, which this code is
 * based on, can be found here:
 * http://www.agner.org/optimize/asmlib.zip
 *
 * This optimization for A/B is helpful if the divisor B is mostly
 * runtime invariant. The reciprocal of B is calculated in the
 * slow-path with reciprocal_value(). The fast-path can then just use
 * a much faster multiplication operation with a variable dividend A
 * to calculate the division A/B.
 */


struct reciprocal_value64 {
	u64 m;
	u8 sh1, sh2;
};

static inline struct reciprocal_value64 reciprocal_value64(u64 d)
{
	struct reciprocal_value64 R;
	int l;

	l = fls64(d - 1);
	
#ifdef WIN32
		uint128 v1;
		v1.Lo = (1ULL << l) - d;v1.Hi=0;
		uint128 v2;
		v2.Hi = 1;v2.Lo = 0;
	
		uint128 v;
		mult128(v1,v2,&v);
		divmod128by64(v.Hi,v.Lo,d,&v.Hi,&v.Lo);
		Increment(&v);	
		R.m = (u64)v.Hi;
#else
    __uint128_t m;
    m = (((__uint128_t)1 << 64) * ((1ULL << l) - d));
    m /= d;
	  ++m;
	  R.m = (u64)m;
#endif
	
	R.sh1 = helpermin(l, 1);
	R.sh2 = helpermax(l - 1, 0);

	return R;
}

static inline u64 reciprocal_divide64(u64 a, struct reciprocal_value64 R)
{
#ifdef WIN32
    uint128 v;
		mult64to128(a,R.m,&v.Hi,&v.Lo);
		u64 t = v.Hi;
#else
    u64 t = (u64)(((__uint128_t)a * R.m) >> 64);
#endif

	return (t + ((a - t) >> R.sh1)) >> R.sh2;
}

static __always_inline uint64_t reciprocal_remainder64(uint64_t A, uint64_t B, struct reciprocal_value64 R)
{
	uint64_t div, mod;

	div = reciprocal_divide64(A, R);
	mod = A - (uint64_t) (div * B);
	if (mod >= B) mod -= B;
	return mod;
}

#endif /* RECIPROCAL_DIV64_H */

