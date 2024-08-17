#ifndef MITM_PROBLEM
#define MITM_PROBLEM

#include "common.hpp"

namespace mitm {
/*
 * Provides a function f : {0, 1}^n -> {0, 1}^n and an optional predicate P.
 *
 * The goal is to find x != y s.t. f(x) == f(y) and P(x, y).
 *
 * The problem instance can contain extra data, e.g. if the goal consists
 * in finding H(prefix || x) == H(prefix || y) with x != y, then the Problem
 * could contain prefix.
 */
class AbstractCollisionProblem {
public:
	int n;           /* size of the domain, in bits */

	u64 f(u64 x) const;
  
	/* assuming that f(x0) == f(x1) and x0 != x1, is (x0, x1) an acceptable outcome? */
	bool is_good_pair(u64 x0, u64 x1) const
	{
		return true;    // by default, yes.
	}
};

/*
 * Provides two functions f, g : {0, 1}^n -> {0, 1}^n and an optional predicate P.
 *
 * The goal is to find x, y s.t. f(x) == g(y) and P(x, y).
 *
 */
class AbstractClawProblem {
public:
	int n;           /* size of the domain, in bits */
  
	u64 f(u64 x) const;
	u64 g(u64 y) const;
  
	/* assuming that f(x0) == g(x1), is (x0, x1) an acceptable outcome? */
	bool is_good_pair(u64 x0, u64 x1) const
	{ 
		return true;    // by default, yes.
	}
};

}
#endif