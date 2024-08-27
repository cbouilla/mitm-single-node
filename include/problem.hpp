#ifndef MITM_PROBLEM
#define MITM_PROBLEM

#include "common.hpp"

namespace mitm {
/*
 * Provides a function f : {0, 1}^n -> {0, 1}^m and an optional predicate P.
 *
 * The goal is to find x != y s.t. f(x) == f(y) and P(x, y).
 *
 * The problem instance can contain extra data, e.g. if the goal consists
 * in finding H(prefix || x) == H(prefix || y) with x != y, then the Problem
 * could contain prefix.
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
  	 * if a vectorized implementation is available, set vlen to the right size 
  	 * and override this function without changing its behavior.
  	 */
	void vf(const u64 x[], u64 y[]) const
	{
		for (int i = 0; i < vlen; i++)
			y[i] = f(x[i]);
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
	int n;           /* size of both domains (input), in bits */
	int m;           /* size of the common range (output), in bits */
	static constexpr int vlen = 1;       /* vector width of the vector implementation */

	/* f, g : {0, 1}^n ---> {0, 1}^m */
	u64 f(u64 x) const { return 0; };
	u64 g(u64 y) const { return 0; };

	/* assuming that f(x0) == g(x1), is (x0, x1) an acceptable outcome? */
	bool is_good_pair(u64 x0, u64 x1) const
	{ 
		return true;    // by default, yes.
	}

  	/* 
  	 * if a vectorized implementation is available, set vlen to the right size 
  	 * and override this functions without changing its behavior.
  	 */
	void vfg(const u64 x[], u64 y[], u64 z[]) const
	{
		for (int i = 0; i < vlen; i++) {
			y[i] = f(x[i]);
			z[i] = g(x[i]);
		}
	}
};
}
#endif