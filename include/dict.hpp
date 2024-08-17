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
	struct __attribute__ ((packed)) pcs_entry {
		u32 key;            // consider that it could be avoided
		u64 start;
		u32 length;         // consider that it could be avoided
	};

	const u64 n_slots;     /* How many slots a dictionary have */

	vector<struct pcs_entry> A;
  
	static u64 get_nslots(u64 nbytes)
	{
		return nbytes / (sizeof(struct pcs_entry));
	}

	PcsDict(u64 n_bytes) : n_slots(get_nslots(n_bytes))
	{
		A.resize(n_slots);
		flush();
	}

	/*
	 * Reset keys, and counters.
	 */
	void flush()
	{
		for (u64 i = 0; i < n_slots; i++)
			A[i].key = 0xffffffff;;
	}
  
	bool is_empty(const struct pcs_entry &e) const
	{
		return e.length == 0xffffffff;
	}

	// murmur64
	u64 hash(u64 h) const
	{
		h ^= h >> 33;
		h *= 0xff51afd7ed558ccdull;
		h ^= h >> 33;
		h *= 0xc4ceb9fe1a85ec53ull;
		h ^= h >> 33;
		return h;
	}

  	// return (start', length'), maybe
	optional<pair<u64, u64>> pop_insert(u64 end, u64 start, u64 len)
	{
		u64 h = hash(end);
		u64 idx = h % n_slots;

		struct pcs_entry &e = A[idx];

		u32 key = h >> 32;
		if (e.key != key || is_empty(e)) [[likely]] {
			e.key = key;
			e.start = start;
			e.length = (u32) len;
			return std::nullopt;
		}

		u64 prev_start = e.start;
		u64 prev_len = e.length;
		e.start = start;
		e.length = len;
		return optional(pair(prev_start, prev_len));
	}
};

}
#endif //MITM_SEQUENTIAL_DICT_HPP
