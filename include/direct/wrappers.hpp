#ifndef MITM_DIRECT_WRAPPERS
#define MITM_DIRECT_WRAPPERS

#include <cassert>
#include <type_traits>

#include "tools.hpp"
#include "problem.hpp"
#include "direct/params.hpp"

/*
 * The two problem wrappers the direct engine runs against: a claw problem evaluates f in FILL and g in
 * PROBE, a collision problem evaluates f in both; each knows how to verify and accept a match.
 */

namespace mitm::direct {

/******************************** the problem wrappers ********************************/

/*
 * A claw problem, f, g : {0,1}^n -> {0,1}^m.  FILL evaluates f, PROBE evaluates g; a match (x, y) is
 * verified with f(x) and tested with is_good_pair(x, y).
 */
template <class Problem>
class ClawWrapper {
public:
	const Problem &pb;              /* the problem itself */
	const int n;                    /* domain bits */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low n bits */
	bool choice[2][Problem::vlen];  /* what vfg selects, one vector per phase: all f in FILL, all g in PROBE */
	static constexpr int vlen = Problem::vlen;
	static constexpr const char *kind = "claw";    /* what the banner calls this search */
	static constexpr const char *funcs = "f, g";   /* the functions it evaluates */

	ClawWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n))
	{
		static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
		              "problem not derived from mitm::AbstractClawProblem");
		assert(n <= 63 && m <= 64);
		for (int k = 0; k < vlen; k++) {
			choice[FILL][k] = true;
			choice[PROBE][k] = false;
		}
	}

	/* vlen evaluations of the phase's function: f while filling, g while probing */
	void veval(int phase, const u64 x[], u64 y[]) const
	{
		if constexpr (vlen == 1)
			y[0] = (phase == FILL) ? pb.f(x[0]) : pb.g(x[0]);
		else
			pb.vfg(x, choice[phase], y);
	}

	/* is the collision f(x) == g(y) the one we want? */
	bool good(u64 x, u64 y) const
	{
		return pb.is_good_pair(x, y);
	}

	/* one value every rank must agree on: the functions have no uninitialised state */
	u64 self_test(u64 a, u64 b) const
	{
		return pb.f(a & in_mask) ^ pb.g(b & in_mask);
	}
};

/*
 * A collision problem, f : {0,1}^n -> {0,1}^m: both phases evaluate f.  A match is a pair of distinct
 * preimages, which is_good_pair -- symmetric, by contract -- accepts or not.
 */
template <class Problem>
class CollisionWrapper {
public:
	const Problem &pb;              /* the problem itself */
	const int n;                    /* domain bits */
	const int m;                    /* range bits */
	const u64 in_mask;              /* the low n bits */
	static constexpr int vlen = Problem::vlen;
	static constexpr const char *kind = "collision";   /* what the banner calls this search */
	static constexpr const char *funcs = "f";          /* the function it evaluates */

	CollisionWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n))
	{
		static_assert(std::is_base_of<AbstractCollisionProblem, Problem>::value,
		              "problem not derived from mitm::AbstractCollisionProblem");
		assert(n <= 63 && m <= 64);
	}

	void veval(int, const u64 x[], u64 y[]) const
	{
		if constexpr (vlen == 1)
			y[0] = pb.f(x[0]);
		else
			pb.vf(x, y);
	}

	/* is the collision f(x) == f(y) the one we want?  is_good_pair is symmetric, so either order will do */
	bool good(u64 x, u64 y) const
	{
		return x != y && pb.is_good_pair(x, y);
	}

	u64 self_test(u64 a, u64 b) const
	{
		return pb.f(a & in_mask) ^ pb.f(b & in_mask);
	}
};

}
#endif
