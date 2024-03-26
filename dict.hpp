//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_SEQUENTIAL_DICT_HPP
#define MITM_SEQUENTIAL_DICT_HPP
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <sys/types.h>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <map>

namespace mitm {
template <typename K, typename V, typename Problem> 
struct Dict {


  const size_t n_slots; /* How many slots a dictionary have */
  size_t n_elements = 0; /* Number of the actual elements that are in the dict */
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
    std::cout << "Done resizing keys\n";
    
    values.resize(n_slots);
    std::cout << "Done resizing values\n";
    chain_lengths.resize(n_slots);
    std::cout << "Done resizing chain\n";
    /* fill values with zeros */
    std::fill(keys.begin(), keys.end(), 0);
    std::fill(chain_lengths.begin(), chain_lengths.end(), 0);
  }

  /*
   * Reset keys, and counters.
   */
  void flush()
  {
    std::fill(keys.begin(), keys.end(), 0);
    std::fill(chain_lengths.begin(), chain_lengths.end(), 0);
    n_elements = 0;
  }
  
  bool pop_insert(const K& key, /* output of f or g */
		  const V& value, /* input */
		  uint64_t const chain_length0,
		  V& out, /* popped input element */
		  uint64_t& chain_length1,
		  Problem& Pb)
  {
    
    /// Add (key, some bits of value) to dictionary. If the pair removes
    /// another pair from the dictionary (Because they have the same index)
    /// write the removed key to out and return true.

    uint64_t idx = key % n_slots;
    bool flag = false;
    
    /* Found an empty slot, thus we're adding a new element  */
    if (keys[idx] == 0)
      ++n_elements; 
    
    if (keys[idx] == key) [[unlikely]]
      flag = true; /* found a collision */

    /* RECALL: */
    /* <value:=input, key:=output>  in this order since value is usually larger */
    /* popped elements */
    Pb.C.copy(values[idx], out); // out = values[idx];
    Pb.C.copy(value, values[idx]); // values[idx] = value;
    chain_length1 = chain_lengths[idx];
    /* write the new entries */
    keys[idx] = key;
    chain_lengths[idx] = chain_length0;

    return flag; /* flag  == 1, if we pop a pair and its value match with the key value */
  }
};


// /* A wrapper over std::unordered_map used only for debugging */
// template <typename K, typename V, typename Problem> 
// struct Dict_unorderd_map {


//   const size_t n_slots; /* How many slots a dictionary have */
//   size_t n_elements = 0; /* Number of the actual elements that are in the dict */
//   size_t n_bytes;

//   /* <value:=input, key:=output>  in this order since value is usually larger */
//   /* uint64_t for the chain length */
//   std::unordered_map<K, std::pair<V, uint64_t> > dict;

  
//   Dict_unorderd_map(size_t n_bytes) :
//     n_slots( n_bytes / ( sizeof(K) + sizeof(V) + sizeof(uint64_t)) )
//   {}

//   /*
//    * Reset keys, and counters.
//    */
//   void flush()
//   {
//     dict.clear();
//     n_elements = 0;
//   }
  
//   bool pop_insert(const K& key, /* output of f or g */
// 		  const V& value, /* input */
// 		  uint64_t const chain_length0,
// 		  V& out, /* popped input element */
// 		  uint64_t& chain_length1,
// 		  Problem& Pb)
//   {
    

//     /// Add (key, some bits of value) to dictionary. If the pair removes
//     /// another pair from the dictionary (Because they have the same index)
//     /// write the removed key to out and return true.

//     // uint64_t idx = key % n_slots;
//     bool flag = false;
    
//     /* Found an empty slot, thus we're adding a new element  */
//     if (dict[key] == 0)
//       ++n_elements; 
    
//     if (keys[idx] == key) [[unlikely]]
//       flag = true; /* found a collision */


//     /* RECALL: */
//     /* <value:=input, key:=output>  in this order since value is usually larger */
//     /* popped elements */
//     Pb.C.copy(values[idx], out); // out = values[idx];
//     Pb.C.copy(value, values[idx]); // values[idx] = value;
//     chain_length1 = chain_lengths[idx];
//     /* write the new entries */
//     keys[idx] = key;
//     chain_lengths[idx] = chain_length0;

//     return flag; /* flag  == 1, if we pop a pair and its value match with the key value */
//   }
// };
}    

#endif //MITM_SEQUENTIAL_DICT_HPP
