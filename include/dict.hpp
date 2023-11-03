//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_SEQUENTIAL_DICT_HPP
#define MITM_SEQUENTIAL_DICT_HPP
#include <cstddef>
#include <cstdint>
#include <vector>

template <typename C_t, typename Val_t> /* A trick allows us to use the output type as an input type! */
// add more details before you forget and cry later.
struct Dict {
    const size_t n_slots; /* How many slots a dictionary have */
    size_t n_elmenents = 0; /* Number of the actual elements that are in the dict */
    using dict_pair = typename std::pair<C_t, Val_t >; /* <input, value> */
    std::vector<dict_pair> table;

    Dict(size_t nSlots) : n_slots(nSlots){
        dict_pair zero_pair; /* we don't care about the value of the first element */
        zero_pair.second = 0;

        table.resize(nSlots); /* Resize the dictionary */
        /* fill values with zeros */
        std::fill(table.begin(), table.end(), zero_pair);
    }

    auto pop_insert(const C_t& key, /* this is a trick for stayin in one type */
                    const C_t& value, /* in our case the distinguished point */
                    C_t& out, /* this will be edited when an element gets popped */
                    uint64_t (*extract_k_bits)(C_t&, uint64_t k))
                    -> bool
    {
        /// Add (key, some bits of value) to dictionary. If the pair removes
        /// another pair from the dictionary (Because they have the same index)
        /// write the removed key to out and return true.

        uint32_t idx = extract_k_bits(key, sizeof(Val_t) ) % n_slots;
        bool flag = 0;
        ++n_elmenents; /* since we are going to add an element */
        /* probe the dictionary at location idx */
        Val_t val = extract_k_bits(value, sizeof(Val_t)); // todo remove hardcoding
        if (table[idx].second != 0){ /* not an empyt slot */
            out = table[idx].first;
            flag = ( table[idx].second == value);
            --n_elmenents; /* we're kicking an element from the dictionary */
        }
        table[idx].first = key;
        table[idx].second = val;

        return flag; /* flag  == 1, if we pop a pair and its value match with the key value */
    }
};

#endif //MITM_SEQUENTIAL_DICT_HPP
