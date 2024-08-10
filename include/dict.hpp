//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_SEQUENTIAL_DICT_HPP
#define MITM_SEQUENTIAL_DICT_HPP

#include <cstdio>
#include <vector>
#include <optional>

#include "common.hpp"

namespace mitm {

template <typename V> 
struct Dict {
  const u64 n_slots; /* How many slots a dictionary have */
  u64 n_elements = 0; /* Number of the actual elements that are in the dict */
  u64 n_bytes;

  /* <value:=input, key:=output>  in this order since value is usually larger */
  u64 *keys;
  V *values;
  
  static u64 get_nslots(u64 nbytes)
  {
    return nbytes / (sizeof(u64) + sizeof(V));
  }

  Dict(u64 n_bytes) : n_slots(get_nslots(n_bytes))
  {
    keys = (u64 *) malloc(n_slots * sizeof(u64));    
    values = (V *) malloc(n_slots * sizeof(V));
    flush();
  }

  ~Dict()
  {
    free(keys);
    free(values);
  }

  /*
   * Reset keys, and counters.
   */
  void flush()
  {
    #pragma omp parallel for
    for (u64 i = 0; i < n_slots; i++)
      keys[i] = 0;
    n_elements = 0;
  }
  
  std::optional<V> pop_insert(u64 key, const V& value)
  {
    /// Add (key, some bits of value) to dictionary. If the pair removes
    /// another pair from the dictionary (Because they have the same index)
    /// write the removed key to out and return true.
    u64 idx = key % n_slots;

    /* Found an empty slot, thus we're adding a new element  */
    if (keys[idx] == 0) {
        #pragma omp atomic
        n_elements += 1;
    }

    if (keys[idx] != key) [[likely]] {
        values[idx] = value;
        keys[idx] = key;
        return std::optional<V>();
    }

    std::optional<V> result(values[idx]);
    values[idx] = value;
    return result;
  }
};

}
#endif //MITM_SEQUENTIAL_DICT_HPP
