//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_SEQUENTIAL_DICT_HPP
#define MITM_SEQUENTIAL_DICT_HPP
#include <cstddef>
#include <cstdint>
#include <sys/types.h>
#include <vector>


namespace mitm {
template <typename K, typename V> 
struct Dict {


  const size_t n_slots; /* How many slots a dictionary have */
  size_t n_elmenents = 0; /* Number of the actual elements that are in the dict */
  size_t n_bytes;

  /* <value:=input, key:=output>  in this order since value is usually larger */
  std::vector<K> keys;
  std::vector<V> values;
  std::vector<uint64_t> chain_lengths;
  
  Dict(size_t n_bytes) :
    n_slots( n_bytes / ( sizeof(K) + sizeof(V) + sizeof(uint64_t)) )
  {
    /*
     * Create a dictionary that takes approximately n_bytes
     */


    /* we don't care about the value of the first element */


    keys.resize(n_slots);
    values.resize(n_slots);
    chain_lengths.resize(n_slots);
    
    /* fill values with zeros */
    std::fill(keys.begin(), keys.end(), 0);
    std::fill(chain_lengths.begin(), chain_lengths.end(), 0);
  }

  bool pop_insert(const K& key, /* output of f or g */
		  const V& value, /* input */
		  uint64_t const chain_length0,
		  V& out, /* popped input element */
		  uint64_t& chain_length1)
  {
    

    /// Add (key, some bits of value) to dictionary. If the pair removes
    /// another pair from the dictionary (Because they have the same index)
    /// write the removed key to out and return true.

    uint64_t idx = key % n_slots;
    bool flag = false;

    /* Found an empty slot, thus we're adding a new element  */
    if (keys[idx] == 0)
      ++n_elmenents; 
    
    if (keys[idx] == key) [[unlikely]]
      flag = true; /* found a collision */

    /* RECALL: */
    /* <value:=input, key:=output>  in this order since value is usually larger */
    /* popped elements */
    out = values[idx];
    chain_length1 = chain_lengths[idx];
    /* write the new entries */
    values[idx] = value;
    keys[idx] = key;

    chain_lengths[idx] = chain_length0;

    return flag; /* flag  == 1, if we pop a pair and its value match with the key value */
  }
};
  
}    

#endif //MITM_SEQUENTIAL_DICT_HPP
