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