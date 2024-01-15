#ifndef MITM_SEQUENTIAL
#define MITM_SEQUENTIAL
#include <array>
#include <cstddef>
#include <cstdint>
#include<cstring>
#include <ios>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <assert.h>
#include <tuple>
#include <utility>
#include <vector>
#include <algorithm>
#include <omp.h>
#include "include/dict.hpp"
#include "include/memory.hpp"
#include <string>

#include <fstream>
#include <chrono>


namespace mitm
{
using u8  = uint8_t ;
using u16 = uint16_t;
using u32 = uint32_t;
using u64 = uint64_t;
using i8  = int8_t ;
using i16 = int16_t;
using i32 = int32_t;
using i64 = int64_t;


/******************************************************************************/
// independent code
/******************************************************************************/
/* Setting up the problem                                                     */
/******************************************************************************/
/* source : https://gist.github.com/mortenpi/9745042 */

template<class T>
T read_urandom()
{
    union {
        T value;
        char cs[sizeof(T)];
    } u;

    std::ifstream rfin("/dev/urandom");
    rfin.read(u.cs, sizeof(u.cs));
    rfin.close();

    return u.value;
}




inline auto wtime() -> double /* with inline it doesn't violate one definition rule */
{

  auto clock = std::chrono::high_resolution_clock::now();
  auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(clock.time_since_epoch()).count();
  double seconds = nanoseconds / (static_cast<double>(1000000000.0));


  return seconds;

}

 /*
  * Print how long does it took to do n iterations.
  */
inline static void print_interval_time(size_t n)
{
  static double previous_time = -1;
  double current_time = wtime();
  double elapsed = current_time - previous_time;
  
  if (previous_time == -1){
    previous_time = wtime();
    return; /* 1st measurement */
  }


  printf("\r%lu=2^%0.2f iter took:"
	 " %0.2f sec, i.e. %0.2f = 2^%0.2f iter/sec",
	 n, std::log2(n),
	 elapsed, n/elapsed, std::log2(n/elapsed) );
  fflush(stdout);
  previous_time = wtime();
  
  
}
  

template <size_t n>
auto is_equal(std::array<u8, n> arr1,
	      std::array<u8, n> arr2)
  -> bool
{
  for (size_t i = 0; i<n; ++i){
    if (arr1[i] != arr2[i])
      return false;
  }
  return true;
}




/******************************************************************************/

/******************************************************************************/
/* Document for standard implementation                                       */
/******************************************************************************/


/*
 * Generic interface for a PRNG. The sequence of pseudo-random numbers
 * depends on both seed and seq
 */
class PRNG {
  /* todo make it support length */
  /* source : https://gist.github.com/mortenpi/9745042 */
  union {
    u64 value;
    char cs[sizeof(u64)];
  } u;
  
public:

  PRNG() { };
  
  u64 rand(){
    std::ifstream rfin("/dev/urandom");
    rfin.read(u.cs, sizeof(u.cs));
    rfin.close();

    return u.value;
  };
  
};


/*
 * A "domain" extends some type to provide the extra functions we need.
 * An instance of the domain can contain extra information
 * E.g. for small integers mod p, "repr" could be u64 while the domain actually contains p
 * E.g. for points on an elliptic curve, "repr" could be a pair of integers mod p, while the
 *      domain would contain the equation of the curve, etc.
 */
template<class repr>           /* repr must support comparisons, and assignment */
class AbstractDomain {
public:
  static int length;  /* nbytes needed to encode an element */
  static size_t n_elements; /* how many elements in the domain */
  using t = repr;            /* t is the machine representation of elements of the domain */

  template<typename PRNG>
  static void randomize(t &x, PRNG &p);           /* set x to a random value */

  static void randomize(t &x);  /* set x to a random value */
  auto is_equal(t& x, t& y) const -> bool;


  /* get the next element after x.*/
  /* What matters is getting a different element each time, not the order. */
  inline t next(t& x)  const;
  inline void serialize(const t &x, u8 *out) const;   /* write this to out */
  inline void unserialize(const u8 *in, t &x) const;        /* read this from in */
  inline void copy(const t& inp, t& out) const; /* deepcopy inp to out */
  inline int extract_1_bit(const t& inp) const;

  u64 hash(const t &x) const;               /* return some bits from this */
  u64 hash_extra(const t &x) const ;         /* return more bits from this */
};


/*
 * Provides the description of the type A and a function f : A -> A.
 *
 * The problem instance can contain extra data.
 * E.g. In the attack on double-encryption, where the goal is to find
 *      x, y s.t. f(x, a) == g(y, b), the problem should contain (a, b).
 */

template<typename Domain_A, typename Domain_B, typename Domain_C>
class AbstractProblem {
public:
  /* these lines have to be retyped again */
  using C_t = typename Domain_C::t;
  using A_t = typename Domain_A::t;
  using B_t = typename Domain_B::t;

  static const int f_eq_g;
  
  AbstractProblem() {
    // enforce that A is a subclass of AbstractDomain
    static_assert(std::is_base_of<AbstractDomain<typename Domain_A::t>, Domain_A>::value,
		  "A not derived from AbstractDomain");
  }
  
  inline void f(const A_t &x, C_t &y) const;  /* y <--- f(x) */
  inline void g(const B_t &x, C_t &y) const;  /* y <--- g(x) */
  inline void send_C_to_A(C_t& inp_C, A_t& out_A) const;
  inline void send_C_to_B(C_t& inp_C, B_t& out_B) const;

  /* changes the behavior of the two above functions */
  void update_embedding(PRNG& rng); 
};

/******************************************************************************/
/* a user does not need to look at the code below                             */
/******************************************************************************/
template<typename C_t>
inline void swap_pointers(C_t*& pt1,
                          C_t*& pt2){
    /// pt1 will point to what pt2 was pointing at, and vice versa.
    C_t* tmp_pt = pt1;
    pt1 = pt2;
    pt2 = tmp_pt;
}

template <typename Problem>
auto is_serialize_inverse_of_unserialize(Problem Pb, PRNG& prng) -> bool
{
  /// Test that unserialize(serialize(r)) == r for a randomly chosen r
  using C_t = typename Problem::C_t;
  const size_t length = Pb.C.length;

  C_t orig{};
  C_t copy{};
  u8 serial[length];

  size_t n_elements = Pb.C.n_elements;

  const size_t n_tests = std::min(n_elements, static_cast<size_t>(1024));


  for(size_t i = 0; i < n_tests; ++i){
    /* */
    Pb.C.randomize(orig, prng);
    Pb.C.serialize(orig, serial);
    Pb.C.unserialize(serial, copy);

    if (not Pb.C.is_equal(copy, orig)){
      /* todo throw exception */
      std::cout << "Error at testing serialization:\n"
		<< "orig = " << orig << "\n"
		<< "copy = " << copy << "\n";
      return false; /* there is a bug in the adaptationof the code  */
    }
      
  }

  return true;

}


   
/* todo pass the function number. e.g. first version of (f, g) second version */
template <typename Problem>
void iterate_once(typename Problem::C_t& inp,
		  typename Problem::C_t& out,
		  Problem& Pb)
{
  /*
   * Do 1 iteration inp =(f/g)=> out, write the output in the address pointed
   * by out_pt.
   */
  typename Problem::A_t inp_A{};
  typename Problem::B_t inp_B{};

  int f_or_g = Pb.C.extract_1_bit(inp);
  if (f_or_g == 1){
    Pb.send_C_to_A(inp, inp_A);
    Pb.f(inp_A, out);
  }
  else { /* f_or_g == 0 */
    Pb.send_C_to_B(inp, inp_B);
    Pb.g(inp_B, out);
  }
}




inline bool is_distinguished_point(u64 digest, u64 mask)
{  return (0 == (mask & digest) ); }

template<typename Problem>
bool generate_dist_point(typename Problem::C_t& inp0, /* don't change the pointer */
			 typename Problem::C_t*& tmp_inp_pt,
			 typename Problem::C_t*& out_pt,
                         u8* const output_bytes, /* size Problem::Dom_C::length */
			 u64& chain_length, /* write chain lenght here */
			 const i64 difficulty, /* #bits are zero at the beginning */
			 Problem& Pb)
{
  /*
   * Given an input, iterate functions either F or G until a distinguished point
   * is found, save the distinguished point in out_pt and output_bytes
   * Then return `true`. If the iterations limit is passed, returns `false`.
   */


  /* copy the input to tmp, then never touch the inp again! */
  Pb.C.copy(inp0, *tmp_inp_pt);

  const u64 mask = (1LL<<difficulty) - 1;
  u64 digest = 0;
  bool found_distinguished = false;
  

  /* The probability of not finding a distinguished point during this loop is*/
  /* 2^{-some big k} */
  for (i64 i = 0; i < 40*(1LL<<difficulty); ++i){
    iterate_once(*tmp_inp_pt, *out_pt, Pb);
    ++chain_length;

    /* we may get a dist point here */
    digest = Pb.C.hash(*out_pt);
    found_distinguished = is_distinguished_point(digest, mask);
      

    /* unlikely with high values of difficulty */
    if (found_distinguished) [[unlikely]]{
      Pb.C.serialize(*out_pt, output_bytes);
      return true; /* exit the whole function */
    }

    /* swap inp and out */
    swap_pointers(tmp_inp_pt, out_pt);
  }
  return false; /* no distinguished point were found */
}





/* Given two inputs that lead to the same distinguished point,
 * find the earliest collision in the sequence before the distinguished point
 * add a drawing to illustrate this.
 */
template<typename Problem>
bool walk(typename Problem::C_t*& inp0_pt,
	  typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
	  u64 inp0_chain_len,
          typename Problem::C_t*& inp1_pt,
	  typename Problem::C_t*& out1_pt, /* inp1 calculation buffer */
	  u64 inp1_chain_len,
	  Problem& Pb)
{
  /****************************************************************************+
   *            walk the longest sequence until they are equal                 |
   * Two chains that leads to the same distinguished point but not necessarily |
   * have the same length. e.g.                                                |
   *                                                                           |
   * chain1: ----------------x-------o                                         |
   *                        /                                                  |
   *          chain2: ------                                                   |
   *                                                                           |
   * o: is a distinguished point                                               |
   * x: the collision we're looking for                                        |
   *                                                                           |   
   ****************************************************************************/

  /* Both sequences need at least `len` steps to reach disitinguish point. */
  size_t const len = std::min(inp0_chain_len, inp1_chain_len);

  /* move the longest sequence until the remaining number of steps is equal */
  /* to the shortest sequence. */
  for (; inp0_chain_len > inp1_chain_len; --inp0_chain_len){
    iterate_once(*inp0_pt, *out0_pt, Pb);
    swap_pointers(inp0_pt, out0_pt);
  }
  
  for (; inp0_chain_len < inp1_chain_len; --inp1_chain_len){
    iterate_once(*inp1_pt, *out1_pt, Pb);
    swap_pointers(inp1_pt, out1_pt);
  }


  /*****************************************************************************/
  /* now both inputs have equal amount of steps to reach a distinguished point */
  /* both sequences needs exactly `len` steps to reach distinguished point.    */
  for (size_t i = 0; i < len; ++i){
    /* walk them together and check each time if their output are equal     */
    /* return as soon equality is found. The equality could be a robinhood. */
    iterate_once(*inp0_pt, *out0_pt, Pb);
    iterate_once(*inp1_pt, *out1_pt, Pb);
    
    if(Pb.C.is_equal( *out0_pt, *out1_pt ))
      return true; /* They are equal */

    swap_pointers(inp0_pt, out0_pt);
    swap_pointers(inp1_pt, out1_pt);
  }
  return false;
}




/*
 * Inputs normall live in A or B, however, the collision only uses C.
 * This function sends the two input to A AND B, if it's not possible, e.g. we can only return 
 */
template <typename Problem>
bool send_2_A_and_B(typename Problem::C_t& inp0_C,
		    typename Problem::C_t& inp1_C,
		    typename Problem::A_t& inp_A,
		    typename Problem::B_t& inp_B,
		    Problem& Pb)
{
  /* when f is same as g, then there is no point in distinguishingg between */
  /* the two functions. */
  if constexpr(Pb.f_eq_g == 1){
    Pb.send_C_to_A(inp0_C, inp_A);
    Pb.send_C_to_B(inp1_C, inp_B);
    return true;
  }
  
  /* otherwise, the collision has to be between f and g */
  int f_or_g = Pb.C.extract_1_bit(inp0_C);
  if (f_or_g == 1){ 
    if (Pb.C.extract_1_bit(inp1_C) == 0){ /* a sensible case*/
      Pb.send_C_to_A(inp0_C, inp_A);
      Pb.send_C_to_B(inp1_C, inp_B);
    }
  } else { /* Pb.C.extract_1_bit(inp0_C) == 0*/
    if (Pb.C.extract_1_bit(inp1_C) == 1){ /* a sensible case*/
      // good case
      Pb.send_C_to_A(inp1_C, inp_A);
      Pb.send_C_to_B(inp0_C, inp_B);
      return true;
    }
  }
  return false;
}


/*
 * return false if the two inputs leads to a robinhood or collision on the same
 * function (when f =/= g). todo: walk should also return false when it could
 * not get a collisison
 */
template <typename Problem >
bool treat_collision(typename Problem::C_t*& inp0_pt,
		     typename Problem::C_t*& out0_pt, /* inp0 calculation buffer */
		     const u64 inp0_chain_len,
		     typename Problem::C_t*& inp1_pt,
		     typename Problem::C_t*& out1_pt, /* inp1 calculation buffer */
		     const u64 inp1_chain_len,
                     std::vector< std::pair<typename Problem::A_t,
		                            typename Problem::B_t> >& container,
		     Problem& Pb)
{
  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;


  /****************************************************************************+
   *            walk the longest sequence until they are equal                 |
   * Two chains that leads to the same distinguished point but not necessarily |
   * have the same length. e.g.                                                |
   *                                                                           |
   * chain1: ----------------x-------o                                         |
   *                        /                                                  |
   *          chain2: ------                                                   |
   *                                                                           |
   * o: is a distinguished point                                               |
   * x: the collision we're looking for                                        |
   *                                                                           |   
   ****************************************************************************/
  /* walk inp0 and inp1 just before `x` */
  /* i.e. iterate_once(inp0) = iterate_once(inp1) */
  /* return false when walking the two inputs don't collide */
  bool found_collision = walk<Problem>(inp0_pt,
				       out0_pt,
				       inp0_chain_len,
				       inp1_pt,
				       out1_pt,
				       inp1_chain_len,
				       Pb);

  /* The two inputs don't lead to the same output */
  if (not found_collision) 
    return false;

  /* If found a robinhood, don't  */
  if ( Pb.C.is_equal(*inp0_pt, *inp1_pt) )
    return false; 

  /* send inp{0,1} to inp_A and inp_B. If it is not possible, return false.*/
  A_t inp_A{};
  B_t inp_B{};


  /* when f=/= g the one of the input should an A input while the other is B's */
  /* If this is not satisfied return */
  bool is_potential_collision = send_2_A_and_B(*inp0_pt, *inp1_pt,
					       inp_A,     inp_B,
					       Pb);
  if (not is_potential_collision)
    return false; /* don't add this pair */
  
  std::pair<A_t, B_t> p{inp_A, inp_B};
  container.push_back(std::move(p)); 
  return true;
}


template <typename Problem, typename C_t>
void apply_f(C_t& inp, Problem& Pb){
  C_t out{};

  Pb.f(inp, out);
  std::cout << "f(inp) = " << out << "\n";
}

template <typename Problem, typename C_t>
void apply_g(C_t& inp, Problem& Pb){
  C_t out{};

  Pb.g(inp, out);
  std::cout << "g(inp) = " << out << "\n";
}



/* todo start from here */
template<typename Problem>
auto collision(Problem& Pb) -> std::pair<typename Problem::C_t, typename Problem::C_t>
{
  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;
  PRNG rng_urandom;

  
  /* Sanity Test: */
  std::cout << "unserial(serial(.)) =?= id(.) : "
	    << is_serialize_inverse_of_unserialize<Problem>(Pb, rng_urandom)
	    << "\n";

  


  // --------------------------------- INIT -----------------------------------/
  size_t n_bytes = 0.5*get_available_memory();
  
  std::cout << "Going to use "
	    << std::dec << n_bytes << " bytes = 2^"<< std::log2(n_bytes)
	    << " bytes for dictionary!\n";
  
  Dict<u64, C_t, Problem> dict{n_bytes}; /* create a dictionary */
  std::cout << "Initialized a dict with " << dict.n_slots << " slots = 2^"
	    << std::log2(dict.n_slots) << " slots\n";


  // -----------------------------------------------------------------------------/
  // VARIABLES FOR GENERATING RANDOM DISTINGUISHED POINTS
  int difficulty = 9; // difficulty;
  /* inp/out variables are used as input and output to save one 1 copy */


  /***************************************************************/
  /* when generating a distinguished point we have:              */
  /*  1)   inp0           =f/g=> out0                            */
  /*  2)  (inp1 := out0)  =f/g=> out1                            */
  /*  3)  (inp2 := out1)  =f/g=> out2                            */
  /*            ...                                              */
  /* m+1) (inp_m := out_m) =f/g=> out_m                          */
  /* A distinguished point found at step `m+1`                   */
  /* "We would like to preserve inp0 at the end of calculation." */
  /* In order to save ourselves from copying in each step.       */
  /***************************************************************/

  C_t pre_inp0{}; /* We need to save this value, see above. */

  /* 1st set of buffers: Related to input0 as a starting point */
  /* either tmp0 or  output0 */
  C_t inp0_or_out0_buffer0{};
  C_t inp0_or_out0_buffer1{};

  C_t* inp0_pt    = &inp0_or_out0_buffer0;
  /* Always points to the region that contains the output */
  C_t* out0_pt = &inp0_or_out0_buffer1;
  u64  out0_digest = 0; /* hashed value of the output0 */
  /* Recall: Problem::Dom_C::length = #needed bytes to encode an element of C_t */
  u8 out0_bytes[Pb.C.length];


  /* 2nd set of buffers: Related to input1 as a starting point */
  /* When we potentially find a collision, we need 2 buffers for (inp1, out1) */
  /* We will use the initial value of inp1 once, thus we don't need to gaurd  */
  /* in an another variable untouched */
  C_t inp1_or_out1_buffer0{};
  C_t inp1_or_out1_buffer1{};
  C_t* inp1_pt    = &inp1_or_out1_buffer0;
  /* Always points to the region that contains the output */
  C_t* out1_pt = &inp1_or_out1_buffer1;

  /* Use these variables to print the full collision */
  A_t inp_A{};
  B_t inp_B{};
  u8 inp_A_serial[Pb.A.length];
  u8 inp_B_serial[Pb.B.length];

  /* Store the results of collisions here */
  /* a:A_t -f-> x <-g- b:B_t */ 
  std::vector< std::pair<A_t, B_t> >  collisions_container{};



  /**************************** Collisions counters ***************************/
  /* How many steps does it take to get a distinguished point from  an input */
  size_t chain_length0 = 0;
  size_t chain_length1 = 0;
  
  bool is_collision_found = false;
  size_t n_collisions = 0;
  size_t n_needed_collisions = 1LL<<20;
  
  /* We should have ration 1/3 real collisions and 2/3 false collisions */
  size_t n_robinhoods = 0;
  

  /*********************************************************
   * surprise we're going actually use
   * void iterate_once(typename Problem::C_t* inp_pt,
   *		  typename Problem::C_t* out_pt,
   *		  Iterate_F<Problem>& F,
   *		  Iterate_G<Problem>& G)
   *
   *********************************************************/

  bool found_dist = false;
  size_t n_distinguished_points = 0;
  constexpr size_t interval = (1LL<<15);
  constexpr size_t interval_mod = interval - 1;
  double collision_timer = wtime();
  
  std::cout << "about to enter a while loop\n";
  /*------------------- Generate Distinguished Points ------------------------*/
  // while (n_collisions < 1){
  while (n_collisions < n_needed_collisions){

    /* These simulations show that if 10w distinguished points are generated
     * for each version of the function, and theta = 2.25sqrt(w/n) then ...
     */
    /* update F and G by changing `send_C_to_A` and `send_C_to_B` */    
    Pb.update_embedding(rng_urandom);
    dict.flush();
    /* todo reset dictionary since all values will be useless */

    /* to do change the number of distinguished points before updates */
    /* After generating xy distinguished point change the iteration function */
    for (size_t n_dist_points = 0;
	 n_dist_points < 10*(dict.n_slots);
	 ++n_dist_points){

      is_collision_found = false;
      /* fill the input with a fresh random value. */
      Pb.C.randomize(pre_inp0, rng_urandom);

      chain_length0 = 0; /* DO we need to reset it? */

      found_dist = generate_dist_point<Problem>(pre_inp0,
						inp0_pt,
						out0_pt,
						out0_bytes,
						chain_length0,
						difficulty,
						Pb);
      out0_digest = Pb.C.hash(*out0_pt);
      ++n_distinguished_points;

      if (n_distinguished_points % interval_mod == 0) /* n_dist mod 2^25 =?= 0 */
	print_interval_time(interval);
	
      // std::cout << "After generatign dist point\n"
      // 		<< "Found distinguished point? " << found_dist << "\n"
      // 		<< "inp    = " << pre_inp0 << "\n"
      // 		<< "inp0   = " << *inp0_pt << "\n"
      // 		<< "out0   = " << *out0_pt << "\n"
      // 		<< "digest = 0x" << std::hex <<  out0_digest << "\n"
      // 		<< "chain length = " << std::dec << chain_length0 << "\n"
      // 		<< "-------\n";
      
      if (not found_dist) [[unlikely]]
	continue; /* skip all calculation below and try again  */
      
      ++n_dist_points;

      
      
      is_collision_found = dict.pop_insert(out0_digest, /* key */
					   pre_inp0, /* value  */
					   chain_length0,
					   *inp1_pt, /* save popped element here */
					   chain_length1, /* of popped input */
					   Pb);
      
      if (is_collision_found) [[unlikely]]{
	/* todo show the results of collisions */
	++n_collisions;

        std::cout << "\nA collision is found\n"
		  << "It took " << (wtime() - collision_timer) << " sec\n"
		  << "inp0 (starting point) = " << pre_inp0 << "\n"
		  << "digest0 = 0x" << out0_digest << "\n"
		  << "chain length0 = " << chain_length0 << "\n"
		  << "inp1 (starting point) = " << *inp1_pt << "\n"
	  	  << "chain length1 = " << std::dec << chain_length1 << "\n"
		  << "-------\n";

	collision_timer = wtime();
	/* respect the rule that inp0 doesn't have pointers dancing around it */
	Pb.C.copy(pre_inp0, *inp0_pt); /* (*tmp0_ptO) holds the input value  */

	/* i.e. when f =/= g then one of the inputs has to  correspond to A
	 * and the other has to correspond to B. The order doesn't matter.
	 * todo: it should also neglect robinhood.
	 */
	bool is_potential_coll = treat_collision<Problem>(inp0_pt,
							  out0_pt,
							  chain_length0,
							  inp1_pt,
							  out1_pt,
							  chain_length1,
							  collisions_container,
							  Pb);



	bool real_collision = Pb.C.is_equal(*out0_pt, *out1_pt);
	bool is_robinhood = Pb.C.is_equal(*inp0_pt, *inp1_pt);
	n_robinhoods += is_robinhood;

	if (is_potential_coll){ 
	  std::cout << "After treating collision\n"
		    << "inp0 = " << *inp0_pt << "\n"
		    << "out0 = " << *out0_pt << "\n"
		    << "inp1 = " << *inp1_pt << "\n"
		    << "out1 = " << *out1_pt << "\n"
		    << "out0 == out1? " << real_collision  << "\n"
		    << "diges0 == digest1? "
		    << "robinhood? " << is_robinhood  << "\n"
		    << "#collisions = " << std::dec << n_collisions << "\n"
		    << "#robinhood = "  << std::dec << n_robinhoods << "\n"
		    << "\n";


	  /* Get the complete inputs as they live in A and B */
	  std::cout << "container length " << collisions_container.size() << "\n"
		    << "is a good collisision? " << is_potential_coll << "\n";
	  Pb.A.serialize(collisions_container.back().first, inp_A_serial);
	  Pb.B.serialize(collisions_container.back().second, inp_B_serial);

	  printf("inp_A = {");
	  for(size_t j = 0; j < Pb.A.length; ++j)
	    printf("0x%02x, ", inp_A_serial[j]);
	  puts("};");

	  printf("inp_B = {");
	  for(size_t j = 0; j < Pb.B.length; ++j)
	    printf("%02x, ", inp_B_serial[j]);
	  puts("};\n________________________________________\n");
	}
      }
    }
  }
  /* end of work */
  return std::pair<C_t, C_t>(*inp0_pt, *inp1_pt); // todo wrong values
}


// to use parallel sort sort
// install tbb lib
// sudo apt install libtbb-dev
// compiling
// g++ -ggdb -O0 -std=c++17 demos/speck32_demo.cpp -o speck32_demo -ltbb
// g++ -flto -O3 -std=c++17  -fopenmp demos/speck32_demo.cpp -o speck32_demo -ltbb

  
}

#endif
