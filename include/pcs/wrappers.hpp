#ifndef MITM_PCS_WRAPPERS
#define MITM_PCS_WRAPPERS

#include <cassert>
#include <type_traits>

#include "tools.hpp"
#include "problem.hpp"

/*
 * The problem seen as one random function of the range onto itself, which is what the trails walk: a
 * version i of the mixing turns a collision or a claw problem into x -> f(mix(i, x)).  A round draws a
 * fresh i, so a run is a sequence of independent functions and the pair is found in one of them.
 */

namespace mitm::pcs {

/*
 * A collision problem, f : {0,1}^n -> {0,1}^m with m >= n, as x -> f(mix(i, x)) on {0,1}^m.  mix is
 * affine modulo 2^n, hence a bijection when the range and the domain have the same size; when the range
 * is the larger one it has to compress, and folds the high bits down so that every input bit is used.
 */
template <class Problem>
class CollisionWrapper {
public:
	const Problem &pb;              /* the problem itself */
	const int n;                    /* domain bits */
	const int m;                    /* range bits: the space the trails walk in */
	const u64 in_mask;              /* the low n bits */
	const u64 out_mask;             /* the low m bits */
	const int fold_shift;           /* how far the range's high bits are folded onto the domain */
	const u64 fold_mask;            /* 0 when the two are the same size, and mix is then a bijection */
	static constexpr int vlen = Problem::vlen;
	static constexpr const char *kind = "collision";   /* what the banner calls this search */
	static constexpr const char *funcs = "f";          /* the function it evaluates */

	CollisionWrapper(const Problem &pb)
		: pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m)),
		  fold_shift((pb.n >= 64) ? 63 : pb.n), fold_mask((pb.m > pb.n) ? ~0ull : 0)
	{
		static_assert(std::is_base_of<AbstractCollisionProblem, Problem>::value,
		              "problem not derived from mitm::AbstractCollisionProblem");
		assert(m <= 64);
		assert(n <= m);

		/* self-test: vmixf agrees with mixf */
		PRNG vprng;
		u64 i = vprng.rand() & out_mask;
		u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		for (int j = 0; j < vlen; j++)
			x[j] = vprng.rand() & out_mask;
		vmixf(i, x, y);
		for (int j = 0; j < vlen; j++)
			assert(y[j] == mixf(i, x[j]));
	}

	/* version i of the mixing: the range onto the domain */
	u64 mix(u64 i, u64 x) const
	{
		u64 y = x * (i | 1) + i;
		return (y ^ ((y >> fold_shift) & fold_mask)) & in_mask;
	}

	u64 mixf(u64 i, u64 x) const
	{
		return pb.f(mix(i, x));
	}

	void vmixf(u64 i, const u64 x[], u64 r[]) const
	{
		/* careful: vlen can be more than one SIMD vector */
		u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		for (int j = 0; j < vlen; j++)
			y[j] = mix(i, x[j]);
		/* a scalar problem provides no vf(): call f() straight away */
		if constexpr (vlen == 1)
			r[0] = pb.f(y[0]);
		else
			pb.vf(y, r);
	}

	/*
	 * Is the collision of the mixed function the pair we want?  Two points that the mixing itself sends
	 * to the same preimage collide under f for no reason of f's own, and are no collision at all.
	 */
	bool mix_good_pair(u64 i, u64 a, u64 b) const
	{
		u64 x0 = mix(i, a);
		u64 x1 = mix(i, b);
		if (x0 == x1)
			return false;
		return pb.is_good_pair(x0, x1);
	}

	/* the answer, in the problem's own terms */
	pair<u64, u64> unmix(u64 i, u64 a, u64 b) const
	{
		return pair(mix(i, a), mix(i, b));
	}

	/* one value every rank must agree on: the function has no uninitialised state */
	u64 self_test(u64 a, u64 b) const
	{
		return mixf(a & out_mask, b & out_mask);
	}
};


/****************************************************************************************/

/*
 * A claw problem with |Domain| == |Range|, as one function: the mixed value picks f or g by one of its
 * bits, so that a collision of the mixed function is a claw whenever the two sides differ.
 */
template <class Problem>
class EqualSizeClawWrapper {
public:
	const Problem &pb;              /* the problem itself */
	const int n;                    /* domain bits */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low n bits */
	const u64 out_mask;             /* the low m bits */
	const u64 choice_mask;          /* the bit of the mixed value that picks f or g */
	static constexpr int vlen = Problem::vlen;
	static constexpr const char *kind = "claw";     /* what the banner calls this search */
	static constexpr const char *funcs = "f, g";    /* the functions it evaluates */

	EqualSizeClawWrapper(const Problem &pb)
		: pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m)),
		  choice_mask(1ull << (pb.m - 1))
	{
		static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
		              "problem not derived from mitm::AbstractClawProblem");
		assert(m <= 64);
		assert(pb.n == pb.m);

		/* self-test: vmixf agrees with mixf */
		PRNG vprng;
		u64 i = vprng.rand() & out_mask;
		u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		for (int j = 0; j < vlen; j++)
			x[j] = vprng.rand() & out_mask;
		vmixf(i, x, y);
		for (int j = 0; j < vlen; j++)
			assert(y[j] == mixf(i, x[j]));
	}

	/* pick either f() or g() */
	bool choose(u64 i, u64 x) const
	{
		return (x * (i | 1)) & choice_mask;
	}

	u64 mix(u64 i, u64 x) const
	{
		return i ^ x;
	}

	u64 mixf(u64 i, u64 x) const
	{
		u64 y = mix(i, x);
		if (choose(i, x))
			return pb.f(y);
		else
			return pb.g(y);
	}

	void vmixf(u64 i, const u64 x[], u64 r[]) const
	{
		/* careful: vlen can be more than one SIMD vector */
		u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		bool choices[vlen];
		for (int j = 0; j < vlen; j++) {
			y[j] = mix(i, x[j]);
			choices[j] = choose(i, x[j]);
		}
		/* a scalar problem provides no vfg(): call f() / g() straight away */
		if constexpr (vlen == 1)
			r[0] = choices[0] ? pb.f(y[0]) : pb.g(y[0]);
		else
			pb.vfg(y, choices, r);
	}

	/* the answer, in the problem's own terms: f's side first, as is_good_pair reads it */
	pair<u64, u64> unmix(u64 i, u64 a, u64 b) const
	{
		u64 x0 = choose(i, a) ? a : b;
		u64 x1 = choose(i, b) ? a : b;
		assert(choose(i, x0));
		assert(not choose(i, x1));
		return pair(mix(i, x0), mix(i, x1));
	}

	bool mix_good_pair(u64 i, u64 a, u64 b) const
	{
		if (choose(i, a) == choose(i, b))
			return false;               /* both on the same function: a collision, but no claw */
		auto [x0, x1] = unmix(i, a, b);
		return pb.is_good_pair(x0, x1);
	}

	/* one value every rank must agree on: the functions have no uninitialised state */
	u64 self_test(u64 a, u64 b) const
	{
		return mixf(a & out_mask, b & out_mask);
	}
};


/*
 * A claw problem with |Range| > |Domain|, as one function of the range onto itself: the mixing folds the
 * range's high bits down, keeping one of them to pick f or g.
 */
template <class Problem>
class LargerRangeClawWrapper {
public:
	const Problem &pb;              /* the problem itself */
	const int n;                    /* domain bits: pb.n + 1, the choice bit included */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low pb.n bits */
	const u64 out_mask;             /* the low m bits */
	const u64 choice_mask;          /* the bit of the mixed value that picks f or g */
	static constexpr int vlen = Problem::vlen;
	static constexpr const char *kind = "claw";     /* what the banner calls this search */
	static constexpr const char *funcs = "f, g";    /* the functions it evaluates */

	LargerRangeClawWrapper(const Problem &pb)
		: pb(pb), n(pb.n + 1), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m)),
		  choice_mask(1ull << (pb.n + 1))
	{
		static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
		              "problem not derived from mitm::AbstractClawProblem");
		assert(m <= 64);
		assert(n <= m);

		/* self-test: vmixf agrees with mixf */
		PRNG vprng;
		u64 i = vprng.rand() & out_mask;
		u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		for (int j = 0; j < vlen; j++)
			x[j] = vprng.rand() & out_mask;
		vmixf(i, x, y);
		for (int j = 0; j < vlen; j++)
			assert(y[j] == mixf(i, x[j]));
	}

	u64 full_mix(u64 i, u64 x) const
	{
		u64 y = x * (i | 1);
		return (y ^ (y >> n));
	}

	/* pick either f() or g() */
	bool choose(u64 i, u64 x) const
	{
		return full_mix(i, x) & choice_mask;
	}

	u64 mix(u64 i, u64 x) const
	{
		return full_mix(i, x) & in_mask;
	}

	u64 mixf(u64 i, u64 x) const
	{
		u64 z = full_mix(i, x);
		if (z & choice_mask)
			return pb.f(z & in_mask);
		else
			return pb.g(z & in_mask);
	}

	void vmixf(u64 i, const u64 x[], u64 r[]) const
	{
		/* careful: vlen can be more than one SIMD vector */
		u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
		bool choice[vlen];
		for (int j = 0; j < vlen; j++) {
			y[j] = mix(i, x[j]);
			choice[j] = choose(i, x[j]);
		}
		/* a scalar problem provides no vfg(): call f() / g() straight away */
		if constexpr (vlen == 1)
			r[0] = choice[0] ? pb.f(y[0]) : pb.g(y[0]);
		else
			pb.vfg(y, choice, r);
	}

	/* the answer, in the problem's own terms: f's side first, as is_good_pair reads it */
	pair<u64, u64> unmix(u64 i, u64 a, u64 b) const
	{
		u64 x0 = choose(i, a) ? a : b;
		u64 x1 = choose(i, b) ? a : b;
		assert(choose(i, x0));
		assert(not choose(i, x1));
		return pair(mix(i, x0), mix(i, x1));
	}

	bool mix_good_pair(u64 i, u64 a, u64 b) const
	{
		if (choose(i, a) == choose(i, b))
			return false;               /* both on the same function: a collision, but no claw */
		auto [x0, x1] = unmix(i, a, b);
		return pb.is_good_pair(x0, x1);
	}

	/* one value every rank must agree on: the functions have no uninitialised state */
	u64 self_test(u64 a, u64 b) const
	{
		return mixf(a & out_mask, b & out_mask);
	}
};

}
#endif
