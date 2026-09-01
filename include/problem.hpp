#ifndef MITM_PROBLEM
#define MITM_PROBLEM

#include <cassert>

#include "tools.hpp"

namespace mitm {
/*
 * Provides a function f : {0, 1}^n -> {0, 1}^m and an optional predicate P.
 *
 * The goal is to find x != y s.t. f(x) == f(y) and P(x, y).
 *
 * The problem instance can contain extra data, e.g. if the goal consists
 * in finding H(prefix || x) == H(prefix || y) with x != y, then the Problem
 * could contain prefix.
 *
 * Derive from this and provide f().  There is no default implementation to inherit:
 * these declarations describe the interface the engine expects, and a missing member
 * is a compile or link error rather than a function that quietly returns nothing.
 */
class AbstractCollisionProblem {
public:
	int n;           /* size of the domain (input), in bits */
	int m;           /* size of the range  (output, in bits */
	static constexpr int vlen = 1;       /* vector width of the vector implementation */

	/* f : {0, 1}^n ---> {0, 1}^m */
	u64 f(u64 x) const;

	/* assuming that f(x0) == f(x1) and x0 != x1, is (x0, x1) an acceptable outcome? */
	bool is_good_pair(u64 x0, u64 x1) const
	{
		return true;    // by default, yes.
	}

	/*
	 * ONLY if a vectorized implementation is available: set vlen to its width and
	 * provide this, without changing its behavior.  With vlen == 1 the engine calls
	 * f() directly and this is never needed.
	 */
	void vf(const u64 x[], u64 y[]) const;
};

/*
 * Provides two functions f, g : {0, 1}^n -> {0, 1}^n and an optional predicate P.
 *
 * The goal is to find x, y s.t. f(x) == g(y) and P(x, y).
 *
 */
class AbstractClawProblem {
public:
	int n;           /* size of both domains (input), in bits */
	int m;           /* size of the common range (output), in bits */
	static constexpr int vlen = 1;       /* vector width of the vector implementation */

	/* f, g : {0, 1}^n ---> {0, 1}^m */
	u64 f(u64 x) const;
	u64 g(u64 y) const;

	/* assuming that f(x0) == g(x1), is (x0, x1) an acceptable outcome? */
	bool is_good_pair(u64 x0, u64 x1) const
	{
		return true;    // by default, yes.
	}

	/*
	 * ONLY if a vectorized implementation is available: set vlen to its width and
	 * provide this, without changing its behavior.  With vlen == 1 the engine calls
	 * f() / g() directly and this is never needed.
	 */
	void vfg(const u64 x[], const bool choice[], u64 y[]) const;
};
}
#endif
