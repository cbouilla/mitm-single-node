//
// Created by ahmed on 23/10/23.
//

#ifndef MITM_SEQUENTIAL_DICT_HPP
#define MITM_SEQUENTIAL_DICT_HPP
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <sys/types.h>
#include <vector>
#include <iostream>
#include <random>

namespace mitm {
template <typename K, typename V, typename Problem> 
struct Dict {


  const size_t n_slots; /* How many slots a dictionary have */
  /* How many keys are there? since we may have 1 key with different values */
  const size_t bucket_size = 8; /* How many elements does a bucket hold */
  const size_t n_keys; 
  size_t n_elements = 0; /* Number of the actual elements that are in the dict */
  // size_t n_bytes;
  std::minstd_rand simple_rand;
  
  
  /* <value:=input, key:=output>  in this order since value is usually larger */
  std::vector<K> keys;
  std::vector<V> values;
  std::vector<uint64_t> chain_lengths;
  
  Dict(size_t n_bytes, size_t bucket_size) :
    n_slots( n_bytes / ( sizeof(K) + bucket_size*sizeof(V) + sizeof(uint64_t)) ),
    bucket_size(bucket_size),
    n_keys(n_slots/bucket_size)
  {
    
    /*
     * Create a dictionary that takes approximately n_bytes
     */

    /* we don't care about the value of the first element */
    keys.resize(n_keys);
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

  /* return a random number between 0 and bucket_size - 1 */
  inline int random_index() {   return simple_rand() % bucket_size;  }

  /* Copy `n_elements` consecutive values from the bucket that its beginning
   * is found at `start_idx` except the value that have the index = `hole_index`
   * EXAMPLE:
   * hole_idx = 2, n_elements = 3, imagine the bucket we are looking at is
   * [x, x, o, x, y, y], we will copy all the x elements, o is the hole element,
   * and we don't care about the y values.
   */
  void copy_holed_bucket(std::vector<V>& out_values, /* popped input element */
			 std::vector<uint64_t>& out_chain_lengths,
			 size_t start_idx,
			 size_t hole_idx,
			 size_t n_elements, /* to be copied */
			 Problem& Pb
			  ) const
  {

    /* assert n_elements <= bucket_size */
    assert(n_elements <= bucket_size);
    
    /* copy elements left to the hole */
    for (int j = 0; j < hole_idx; ++j){
      Pb.C.copy(values[j + start_idx], out_values[j] ); // values[idx] = value;
      out_chain_lengths[j] = chain_lengths[j + start_idx];
    }
    /* copy elements right to the hole */
    for (int j = hole_idx + 1; j <= n_elements; ++j) {
      Pb.C.copy(values[j + start_idx], out_values[j] ); // values[idx] = value;
      out_chain_lengths[j] = chain_lengths[j + start_idx];
    }
  }
  
  /* Return how many values were found with the same inquired key */
  int pop_insert(const K& key, /* output of f or g */
		  const V& value, /* input */
		  uint64_t const chain_length,
		  std::vector<V>& out_values, /* popped input element */
		  std::vector<uint64_t>& out_chain_lengths,
		  Problem& Pb)
  {
    
    /// Add (key, some bits of value) to dictionary. If the pair removes
    /// another pair from the dictionary (Because they have the same index)
    /// write the removed key to out and return true.

    /* How many elements in a bucket */
    int n_values_in_bucket = 0;
    
    uint64_t idx = key % n_keys; /* where to look in the dictionary? */
    /* it will change as soon as we found an empty slot in a bucket*/
    bool bucket_was_full = true;
    
    /* Find an empty slot in the bucket. Count how many elements are there */
    for (int j = 0; j < bucket_size; ++j){
      /* Found an empty slot */
      if (chain_lengths[j + idx*bucket_size] == 0){
	bucket_was_full = false;
	chain_lengths[j + idx*bucket_size] = chain_length;
	Pb.C.copy(value, values[j + idx*bucket_size]); // values[idx] = value;
	break; /* no need to continue */
      }
      ++n_values_in_bucket;
    }

    /* where the element was stored in the bucket. If the bucket was full, then
     * the hole is outside the bucket. */
    size_t hole_idx = n_values_in_bucket;

    /* Copy exactly the values that already exist in bucket except the hole */
    copy_holed_bucket(out_values,
		      out_chain_lengths,
		      idx*bucket_size,
		      hole_idx,
		      n_values_in_bucket,
		      Pb);

    /* if the element was not inserted since the bucket was full. Pick a
     * random element in the bucket to kicked out. */
    if (bucket_was_full){
      hole_idx = random_index();
      chain_lengths[hole_idx + idx*bucket_size] = chain_length;
      Pb.C.copy(value, values[hole_idx + idx*bucket_size]); // values[idx] = value;
    }

    /* number of elements we've extracted from the dictionary to:
     * out_values and out_chain_lenths. */
    return n_values_in_bucket;
  }
};


}    

#endif //MITM_SEQUENTIAL_DICT_HPP
