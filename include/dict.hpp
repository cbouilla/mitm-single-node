//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_SEQUENTIAL_DICT_HPP
#define MITM_SEQUENTIAL_DICT_HPP

#include "tools.hpp"

// various dictionnaries

namespace mitm {

/*
 * this is a "classic" hash table for 64-bit key-value pairs, with linear probing.  
 * No false negatives, some false positives.  12 bytes per entry. 
 */
class CompactDict {
public:
    const u64 n_slots;     /* How many slots a dictionary have */
    struct __attribute__ ((packed)) entry { u32 k; u64 v; };

    vector<struct entry> A;

    CompactDict(u64 n_slots) : n_slots(n_slots)
    {
        A.resize(n_slots, {0xffffffff, 0});
    }

    void insert(u64 key, u64 value)
    {
        u64 h = (key ^ (key >> 32)) % n_slots;
        for (;;) {
            if (A[h].k == 0xffffffff)
                break;
            h += 1;
            if (h == n_slots)
                h = 0;
        }
        A[h].k = key % 0xfffffffb;
        A[h].v = value;
    }

    // return possible values matching this key
    // TODO: replace this by a custom iterator
    int probe(u64 bigkey, u64 keys[])
    {
        u32 key = bigkey % 0xfffffffb;
        u64 h = (bigkey ^ (bigkey >> 32)) % n_slots;
        int nkeys = 0;
        for (;;) {
            if (A[h].k == 0xffffffff)
                return nkeys;
            if (A[h].k == key) {
                keys[nkeys] = A[h].v;
                nkeys += 1;
            }
            h += 1;
            if (h == n_slots)
                h = 0;
        }
    }
};

/*
 * This dictionnary, when probed with the distinguished point at the end of a trail,
 * should provide (if any) the start (and the length) of another distinguished point
 * that has the same end.
 */
class PcsDict {
public:
	u64 m;
	u64 start_mask;
	u64 key_mask;
	const u64 n_slots;     /* size of A */
	vector<u64> A;
  
	static u64 get_nslots(u64 nbytes, u64 forced_multiple)
	{
		u64 w = nbytes / (sizeof(u64));
		return (w / forced_multiple) * forced_multiple;
	}

	PcsDict(u64 m, u64 w) : m(m), n_slots(w)
	{
		start_mask = make_mask(m);
		key_mask = (m == 64) ? 0 : 0xffffffffffffffff << m;
		A.resize(n_slots);
		flush();
	}

	/*
	 * Reset keys, and counters.
	 */
	void flush()
	{
		for (u64 i = 0; i < n_slots; i++)
			A[i] = 0;
	}
  
  	// return (start'), maybe
	optional<u64> pop_insert(u64 end, u64 start)
	{
		u64 idx = end % n_slots;
		u64 key = ((end / n_slots) << m) & key_mask;

		u64 e = A[idx];
		u64 ekey = e & key_mask;
		assert((start & start_mask) == start);
		A[idx] = start ^ key;

		if (ekey != key || e == 0)
			return nullopt;
		else
			return optional(e & start_mask);
	}
};

}
#endif //MITM_SEQUENTIAL_DICT_HPP
