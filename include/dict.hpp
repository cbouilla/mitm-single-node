//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_SEQUENTIAL_DICT_HPP
#define MITM_SEQUENTIAL_DICT_HPP
#include <cstddef>
#include <cstdint>
#include <vector>

template <typename C_t> /* A trick allows us to use the output type as an input type! */
// add more details before you forget and cry later.
struct Dict {
    size_t n_slots; /* How many slots a dictionary have */
    size_t n_elmenents = 0; /* Number of the actual elements that are in the dict */
    using dict_pair = typename std::pair<C_t, uint32_t >; /* <input, value> */
    std::vector<dict_pair> table;

    Dict(size_t nSlots) : n_slots(nSlots){
        dict_pair zero_pair; /* we don't care about the value of the first element */
        zero_pair.second = 0;

        table.resize(nSlots); /* Resize the dictionary */
        /* fill values with zeros */
        std::fill(table.begin(), table.end(), zero_pair);
    }
    auto pop_insert(const C_t& inp,
                    uint8_t* inp_serialized,  // maybe we should not add C_t as a variable
                    const C_t& dist_point,
                    uint8_t* dist_point_serialized,
                    C_t& out)
                    -> bool
    {
        /* todo we should only accept serialized input, or maybe not? */

        uint32_t idx = *(reinterpret_cast<uint32_t*> (&inp_serialized[0])); // dear loard, this is ugly!
        bool flag = 0;
        /* probe the dictionary at location idx */
        uint32_t value = *(reinterpret_cast<uint32_t*> (&dist_point_serialized[0])); /* value TO BE stored in the dictionary */
        if (table[idx].second != 0){ /* not an empyt slot */
            out = table[idx].first;

            flag = ( table[idx].second == value);
        }
        table[idx].first = inp;
        table[idx].second = value;

        return flag; /* flag  == 1, if we pop a pair and its value match with the inp value */
    }
};

#endif //MITM_SEQUENTIAL_DICT_HPP
