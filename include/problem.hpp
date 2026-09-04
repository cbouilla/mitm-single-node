#ifndef MITM_PROBLEM
#define MITM_PROBLEM

#include <cassert>

#include "tools.hpp"

namespace mitm {
/*
 * A collision problem: f : {0, 1}^n -> {0, 1}^m and a predicate on the pair.  The goal is x != y with
 * f(x) == f(y) and is_good_pair(x, y).  Derive and provide f(); nothing here is defined, on purpose: a
 * missing member is a compile or link error rather than a function that quietly returns 0.
 */
class AbstractCollisionProblem {
public:
	int n;           /* size of the domain (input), in bits */
	int m;           /* size of the range (output), in bits */
	static constexpr int vlen = 1;       /* vector width of the vector implementation */

	/* f : {0, 1}^n ---> {0, 1}^m */
	u64 f(u64 x) const;

	/* assuming that f(x0) == f(x1) and x0 != x1, is (x0, x1) an acceptable outcome? */
	bool is_good_pair(u64 x0, u64 x1) const;

	/* only with vlen > 1: f() on vlen inputs at once, same results */
	void vf(const u64 x[], u64 y[]) const;
};

/*
 * A claw problem: f, g : {0, 1}^n -> {0, 1}^m and a predicate on the pair.  The goal is x, y with
 * f(x) == g(y) and is_good_pair(x, y).  Same rules as above.
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
	bool is_good_pair(u64 x0, u64 x1) const;

	/* only with vlen > 1: f() where choice[j] is set, g() elsewhere, same results */
	void vfg(const u64 x[], const bool choice[], u64 y[]) const;
};
}
#endif
