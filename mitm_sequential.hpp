#ifndef MITM_SEQUENTIAL
#define MITM_SEQUENTIAL
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ios>
#include <cmath>
#include <iostream>
#include <iterator>
#include <stdio.h>
#include <assert.h>
#include <pstl/glue_execution_defs.h>
#include <tuple>
#include <utility>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <execution>
#include "include/dict.hpp"
#include <string>

size_t n_f_called{0};
size_t n_g_called{0};
size_t n_robinhood{0};

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
#include <fstream>
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



#include <chrono>
inline auto wtime() -> double /* with inline it doesn't violate one definition rule */
{

  auto clock = std::chrono::high_resolution_clock::now();
  auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(clock.time_since_epoch()).count();
  double seconds = nanoseconds / (static_cast<double>(1000000000.0));


  return seconds;

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
  inline void serialize(const t &x, const u8 *out) const;   /* write this to out */
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
  Domain_A A;
  Domain_B B;
  Domain_C C;
  using C_t = typename Domain_C::t;
  using A_t = typename Domain_A::t;
  using B_t = typename Domain_B::t;

  void f(const A_t &x, C_t &y) const;  /* y <--- f(x) */
  void g(const B_t &x, C_t &y) const;  /* y <--- g(x) */


  AbstractProblem() {
    // enforce that A is a subclass of AbstractDomain
    static_assert(std::is_base_of<AbstractDomain<typename Domain_A::t>, Domain_A>::value,
		  "A not derived from AbstractDomain");
  }

  void send_C_to_A(C_t& inp_C, A_t& out_A) const;
  void send_C_to_B(C_t& inp_C, B_t& out_B) const;
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
  using C_t = typename Problem::Dom_C::t;
  const size_t length = Problem::Dom_C::length;

  C_t orig{};
  C_t copy{};
  std::array<u8, length> serial;

  size_t n_elements = Pb.C.n_elements;

  const size_t n_tests = std::min(n_elements, static_cast<size_t>(1024));


  for(size_t i = 0; i < n_tests; ++i){
    /* */
    Pb.C.randomize(orig, prng);
    Pb.C.serialize(orig, serial);
    Pb.C.unserialize(serial, copy);

    if (copy != orig)
      return false; /* there is a bug in the adaptationof the code  */
  }

  return true;

}


   
/* todo pass the function number. e.g. first version of (f, g) second version */
template <typename Problem>
void iterate_once(typename Problem::Dom_C::t& inp,
		  typename Problem::Dom_C::t& out,
		  Problem& Pb)
{
  /*
   * Do 1 iteration inp =(f/g)=> out, write the output in the address pointed
   * by out_pt.
   */
  typename Problem::Dom_A::t inp_A{};
  typename Problem::Dom_A::t inp_B{};

  int f_or_g = Pb.C.extract_1_bit(inp);
  if (f_or_g == 1){
    Pb.send_C_to_A(inp, inp_A);
    Pb.f(inp, out);
    ++n_f_called;
  }
  else { /* f_or_g == 0 */
    Pb.send_C_to_B(inp, inp_B);
    Pb.f(inp_B, out);
    ++n_g_called;
  }
}




inline bool is_distinguished_point(u64 digest, u64 mask)
{  return (0 == (mask & digest) ); }

template<typename Problem>
bool generate_dist_point(typename Problem::Dom_C::t& inp0, /* don't change the pointer */
			 typename Problem::Dom_C::t*& tmp_inp_pt,
			 typename Problem::Dom_C::t*& out_pt,
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
      std::cout << "Found dist, digest = 0x" << std::hex << digest << "\n";
      Pb.C.serialize(*out_pt, output_bytes);
      return true; /* exit the whole function */
    }

    /* swap inp and out */
    swap_pointers(tmp_inp_pt, out_pt);
  }
  return false; /* no distinguished point were found */
}






template<typename Problem>
void walk(typename Problem::Dom_C::t*& inp0_pt,
	  typename Problem::Dom_C::t*& tmp0_pt, /* inp0 calculation buffer */
	  const u64 inp0_chain_len,
          typename Problem::Dom_C::t*& inp1_pt,
	  typename Problem::Dom_C::t*& tmp1_pt, /* inp1 calculation buffer */
	  const u64 inp1_chain_len,
	  Problem& Pb)
{
  /* Given two inputs that lead to the same distinguished point,
   * find the earliest collision in the sequence before the distinguished point
   * add a drawing to illustrate this.
   */
  using C_t = typename Problem::Dom_C::t;

  size_t const diff_len = std::max(inp0_chain_len, inp1_chain_len)
                        - std::min(inp0_chain_len, inp1_chain_len);

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
  
 
  if (inp0_chain_len > inp1_chain_len){
    for (size_t i = 0; i < diff_len; ++i){
      iterate_once(*inp0_pt, *tmp0_pt, Pb);
      swap_pointers(inp0_pt, tmp0_pt);
    }
   }
 
  if (inp0_chain_len < inp1_chain_len){
    for (size_t i = 0; i < diff_len; ++i){
      iterate_once(*inp1_pt, *tmp1_pt, Pb);
      swap_pointers(inp0_pt, tmp0_pt);
    }
  }

  /****************************************************************************/
  /* now both inputs have equal amount of steps to reach a distinguished point */
  /* both sequences needs exactly `len` steps to reach dist point */
  size_t len = std::min(inp0_chain_len, inp1_chain_len);

  /* Check if the two inputs are equal then we have a robin hood */


  if(Pb.C.is_equal( *inp0_pt, *inp0_pt ))
    return; /* Robinhood */

  
  for (size_t i = 0; i < len; ++i){
    /* walk them together and check each time if their output are equal */
    /* return as soon equality is found */
    iterate_once(*inp0_pt, *tmp0_pt, Pb);
    iterate_once(*inp1_pt, *tmp1_pt, Pb);
    
    if(Pb.C.is_equal( *tmp0_pt, *tmp1_pt ))
      return; /* They are equal */
  }
}




template <typename Problem >
bool treat_collision(typename Problem::Dom_C::t*& inp0_pt,
		     typename Problem::Dom_C::t*& tmp0_pt, /* inp0 calculation buffer */
		     const u64 inp0_chain_len,
		     typename Problem::Dom_C::t*& inp1_pt,
		     typename Problem::Dom_C::t*& tmp1_pt, /* inp1 calculation buffer */
		     const u64 inp1_chain_len,
                     std::vector< std::pair<typename Problem::Dom_C::t,
		                            typename Problem::Dom_C::t> >& container,
		     Problem& Pb)
{
  using A_t = typename Problem::Dom_A::t;
  using B_t = typename Problem::Dom_B::t;
  using C_t = typename Problem::Dom_C::t;


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
  walk<Problem>(inp0_pt,
		tmp0_pt,
		inp0_chain_len,
		inp1_pt,
		tmp1_pt, /* inp1 calculation buffer */
		inp1_chain_len,
	        Pb);
					   

  /* assume copying */
  std::pair<C_t, C_t> p{*inp0_pt, *inp1_pt};
  
  /* todo here we should add more tests */
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
auto collision(Problem& Pb) -> std::pair<typename Problem::Dom_C::t, typename Problem::Dom_C::t>
{
  using A_t = typename Problem::Dom_A::t;
  using B_t = typename Problem::Dom_B::t;
  using C_t = typename Problem::Dom_C::t;
  PRNG rng_urandom;

  
  /* Sanity Test: */
  std::cout << "unserial(serial(.)) =?= id(.) : "
	    << is_serialize_inverse_of_unserialize<Problem>(Pb, rng_urandom)
	    << "\n";

  


  // --------------------------------- INIT -----------------------------------/
  size_t n_bytes = 1LL<<34; /* */
  Dict<u64, C_t> dict{n_bytes}; /* create a dictionary */
  std::cout << "Initialized a dict with " << dict.n_slots << " slots\n";


  // -----------------------------------------------------------------------------/
  // VARIABLES FOR GENERATING RANDOM DISTINGUISHED POINTS
  int difficulty = 4; // difficulty;
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

  C_t inp0{}; /* We need to save this value, see above. */

  /* 1st set of buffers: Related to input0 as a starting point */
  /* either tmp0 or  output0 */
  C_t inp0_or_out0_buffer0{};
  C_t inp0_or_out0_buffer1{};

  C_t* tmp0_pt    = &inp0_or_out0_buffer0;
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



  /* Store the results of collisions here */
  /* a:A_t -f-> x <-g- b:B_t */ 
  std::vector< std::pair<A_t, B_t> >  collisions_container{};
  // TODO FATAL if A_t =/= B_t we will have an error because treat collision uses std::pair<C_t, C_t





  /**************************** Collisions counters ***************************/
  /* How many steps does it take to get a distinguished point from  an input */
  size_t chain_length0 = 0;
  size_t chain_length1 = 0;
  
  bool is_collision_found = false;
  size_t n_collisions = 0;
  size_t n_needed_collisions = 1LL<<20;

  /* We should have ration 1/3 real collisions and 2/3 false collisions */
  size_t false_collisions = 0;
  size_t real_collisions = 0;
  size_t n_dist_points = 0;



  /*********************************************************
   * surprise we're going actually use
   * void iterate_once(typename Problem::Dom_C::t* inp_pt,
   *		  typename Problem::Dom_C::t* out_pt,
   *		  Iterate_F<Problem>& F,
   *		  Iterate_G<Problem>& G)
   *
   *********************************************************/

  bool found_dist = false;
  
  /*------------------- Generate Distinguished Points ------------------------*/
  while (n_collisions < n_needed_collisions){



    /* These simulations show that if 10w distinguished points are generated
     * for each version of the function, and theta = 2.25sqrt(w/n) then ...
     */
    /* update F and G by changing `send_C_to_A` and `send_C_to_B` */    
    Pb.update_embedding(rng_urandom);
    /* todo reset dictionary since all values will be useless */

    /* to do change the number of distinguished points before updates */
    /* After generating xy distinguished point change the iteration function */
    for (size_t n_dist_points = 0; n_dist_points < (1LL<<30); ++n_dist_points){
      is_collision_found = false;
      /* fill the input with a fresh random value. */
      Pb.C.randomize(inp0, rng_urandom);

      chain_length0 = 0; /* DO we need to reset it? */
      std::cout << "==================\n"
		<< "At the beginning\n"
		<< "inp = " << inp0 << "\n"
		<< "out = " << *out0_pt << "(garbage from previous calculations!)\n"
		<< "chain length = " << chain_length0 << "\n"
		<< "-------\n";

      found_dist = generate_dist_point<Problem>(inp0,
						tmp0_pt,
						out0_pt,
						out0_bytes,
						chain_length0,
						difficulty,
						Pb);

      if (not found_dist) [[unlikely]]
	continue; /* skip all calculation below and try again  */
      
      ++n_dist_points;
      out0_digest = Pb.C.hash(*out0_pt);
      
      std::cout << "After generatign dist point\n"
		<< "Found distinguished point? " << found_dist << "\n"
		<< "inp    = " << inp0 << "\n"
		<< "tmp0   = " << *tmp0_pt << "\n"
		<< "out0   = " << *out0_pt << "\n"
		<< "digest = 0x" << std::hex <<  out0_digest << "\n"
		<< "chain length = " << chain_length0 << "\n"
		<< "-------\n";
      

      
      is_collision_found = dict.pop_insert(out0_digest, /* key */
					   inp0, /* value  */
					   *inp1_pt, /* save popped element here */
					   chain_length1 /* of popped input */);
      
      if (is_collision_found) [[unlikely]]{
	/* todo show the results of collisions */
	++n_collisions; 
	bool is_false_collision = false;


	std::cout << "A collision is found\n"
		  << "inp0    = " << inp0 << "\n"
		  << "digest0 = " << out0_digest << "\n"
		  << "chain length0 = " << chain_length0 << "\n"
		  << "inp1    = " << *inp1_pt << "\n"
	  	  << "chain length1 = " << chain_length1 << "\n"
		  << "-------\n";

	
	/* respect the rule that inp0 doesn't have pointers dancing around it */
	Pb.C.copy(inp0, *tmp0_pt); /* (*tmp0_pt) holds the input value  */

	

	
	treat_collision<Problem>(tmp0_pt, /* inp0 */
				 out0_pt, /* inp0 scratch buffer  */
				 chain_length0,
				 inp1_pt,
				 out1_pt,
				 chain_length1,
				 collisions_container,
				 Pb);



	std::cout << "After treating collision\n"
		  << "inp0 = " << *tmp0_pt << "\n"
		  << "out0 = " << *out0_pt << "\n"
		  << "inp1 = " << *inp1_pt << "\n"
		  << "out1 = " << *out1_pt << "\n"
		  << "_______\n";

      }
    }
  }
  /* end of work */
  return std::pair<C_t, C_t>(*tmp0_pt, *inp1_pt); // todo wrong values
}
#endif

// to use parallel sort sort
// install tbb lib
// sudo apt install libtbb-dev
// compiling
// g++ -flto -O3 -std=c++17  -fopenmp demos/speck32_demo.cpp -o speck32_demo -ltbb

// thd0 has out = 9ace,f6f8 and inp_out[i] = ce9a,f8f6 and inp = 1,0 and inp_out[i].first = 1,0
//AFTER unserialize thd0           has out = 9ace,f6f8 and inp_out[i] = ce9a,f8f6 and inp = 1,0 and inp_out[i].first = 1,
