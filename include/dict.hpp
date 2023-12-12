//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_SEQUENTIAL_DICT_HPP
#define MITM_SEQUENTIAL_DICT_HPP
#include <cstddef>
#include <cstdint>
#include <vector>

template <typename K, typename V> 

struct Dict {
  const size_t n_slots; /* How many slots a dictionary have */
  size_t n_elmenents = 0; /* Number of the actual elements that are in the dict */

  /* <value:=input, key:=output>  in this order since value is usually larger */
  using dict_pair = typename std::pair<V, K>; 
  std::vector<dict_pair> table;

  Dict(size_t nSlots) : n_slots(nSlots){
    dict_pair zero_pair;
    /* we don't care about the value of the first element */
    zero_pair.second = 0; /* set key to 0 */

    table.resize(nSlots); /* Resize the dictionary */
    /* fill values with zeros */
    std::fill(table.begin(), table.end(), zero_pair);
  }

  bool pop_insert(const K& key,
		  const V& value,
		  V& out, 
		  uint64_t (*extract_k_bits)(const V&, int k))
  
  {
    /// Add (key, some bits of value) to dictionary. If the pair removes
    /// another pair from the dictionary (Because they have the same index)
    /// write the removed key to out and return true.

    uint64_t idx = extract_k_bits(key, sizeof(V )*8) % n_slots;
    bool flag = false;
    ++n_elmenents; /* since we are going to add an element */

    /* probe the dictionary at location idx */
    if (table[idx].second != 0){ /* not an empyt slot */
      out = table[idx].first;
      flag = ( table[idx].second == value);
      if (flag)
          // std::cout << "dict entry  " << table[idx].second
          // << " val = " << val
          // << " idx = " << idx
          // << ", sizeof(Val) = " << sizeof(V) << "\n";

      --n_elmenents; /* we're kicking an element from the dictionary */
    }
    table[idx].first = key;
    table[idx].second = value;

    return flag; /* flag  == 1, if we pop a pair and its value match with the key value */
  }
};

#endif //MITM_SEQUENTIAL_DICT_HPP
