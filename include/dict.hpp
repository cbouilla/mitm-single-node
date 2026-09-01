//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_DICT
#define MITM_DICT

#include <cassert>

#include "tools.hpp"

namespace mitm {

/*
 * This dictionnary, when probed with the distinguished point at the end of a trail,
 * should provide (if any) the start (and the length) of another distinguished point
 * that has the same end.
 */
class PcsDict {
public:
	u64 jbits, lbits;
	u64 jmask, lmask;
	u64 len_mask;
	u64 key_mask;
	const u64 n_slots;     /* size of A */
	
	vector<u64> A;         // A[i][0:jbits] == j.  A[i][jbits:lbits] == len1.  A[lbits:64] == key bits
  	
	static u64 get_nslots(u64 nbytes, u64 forced_multiple)
	{
		u64 w = nbytes / (sizeof(u64));
		return (w / forced_multiple) * forced_multiple;
	}

	PcsDict(u64 jbits, u64 w) : jbits(jbits), n_slots(w)
	{
		assert(jbits <= 56);
		jmask = make_mask(jbits);
		lmask = make_mask(8);
		lbits = jbits + 8;
		key_mask = (lbits == 64) ? 0 : 0xffffffffffffffff << lbits;
		A.resize(n_slots);
	}

	void flush()
	{
		for (u64 i = 0; i < n_slots; i++)
			A[i] = 0;
	}
  
  	// return (start', len'), maybe. Return len' == 0 if len' is unknown (has been truncated)
	optional<pair<u64, u64>> pop_insert(u64 end, u64 start, u64 len0)
	{
		u64 idx = end % n_slots;
		u64 key = (end / n_slots) << lbits;

		// first save the previous point
		u64 e = A[idx];
		u64 ekey = e & key_mask;
		u64 elen = (e >> jbits) & lmask;

		/*
		 * heuristic modification from the original algorithm:
         * we only overwrite an existing point if we have a longer tail.
         * the alternative is to always insert and forget about it. 
         */
		if (e == 0 || len0 >= elen) {
			if (len0 > lmask)      // saturate the length because be don't use many bits for it.
				len0 = lmask;
			A[idx] = start ^ (len0 << jbits) ^ key;
		}

		if (ekey != key || e == 0)
			return nullopt;
		
		if (elen == lmask)
			elen = 0;

		return optional(pair(e & jmask, elen));
	}
};

}
#endif
