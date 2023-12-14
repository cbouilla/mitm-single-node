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
#include <set>
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
#include <chrono>
inline auto wtime() -> double /* with inline it doesn't violate one definition rule */
{

  auto clock = std::chrono::high_resolution_clock::now();
  auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(clock.time_since_epoch()).count();
  double seconds = nanoseconds / (static_cast<double>(1000000000.0));


  return seconds;

}


template <typename T>
void print_array(T arr){
  /* Would not be a better choice if we just to overload << for std::array? */
  printf("0x");
  for (auto& item: arr)
    printf("%x", item);
  printf("\n");
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


template <typename Problem>
auto is_serialize_inverse_of_unserialize() -> bool
{
  /// Test that unserialize(serialize(r)) == r for a randomly chosen r
  using C_t = typename Problem::C::t;
  const size_t length = Problem::C::length;

  C_t orig{};
  C_t copy{};
  std::array<u8, length> serial;

  size_t n_elements = Problem::C::n_elements;

  const size_t n_tests = std::min(n_elements, static_cast<size_t>(1024));


  for(size_t i = 0; i < n_tests; ++i){
    /* */
    Problem::C::randomize(orig);
    Problem::C::serialize(orig, serial);
    Problem::C::unserialize(serial, copy);

    if (copy != orig)
      return false; /* there is a bug in the adaptationof the code  */
  }

  return true;

}


/******************************************************************************/

/******************************************************************************/
/* Document for standard implementation                                       */
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
/*
 * Generic interface for a PRNG. The sequence of pseudo-random numbers
 * depends on both seed and seq
 */
class AbstractPRNG {
public:
  AbstractPRNG(u64 seed, u64 seq);
  u64 rand();
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

  template<class PRNG>
  static void randomize(t &x, PRNG &p);           /* set x to a random value */

  static void randomize(t &x);  /* set x to a random value */
  auto is_equal(t& x, t& y) const -> bool;


  /* get the next element after x. What matters is getting a different element eac time, not the order. */
  inline auto next(t& x) -> t;
  inline static void serialize(const t &x, const u8 *out);   /* write this to out */
  inline static void unserialize(const u8 *in, t &x);        /* read this from in */
  inline static void copy(const t& inp, t& out); /* deepcopy inp to out */
  inline static auto extract_1_bit(const t& inp) -> int;
  inline static auto extract_k_bits(const t& inp, int k) -> u64;

  static auto hash(const t &x) -> u64 ;               /* return some bits from this */
  static auto hash_extra(const t &x) -> u64 ;         /* return more bits from this */
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
  using C =  Domain_C;
  using C_t = typename C::t;

  using A = Domain_A;
  using A_t = typename A::t;


  using B = Domain_B;
  using B_t = typename B::t;

  static void f(const A_t &x, C_t &y);                /* y <--- f(x) */
  static void g(const B_t &x, C_t &y);                /* y <--- g(x) */


  AbstractProblem() {
    // enforce that A is a subclass of AbstractDomain
    static_assert(std::is_base_of<AbstractDomain<typename A::t>, A>::value,
		  "A not derived from AbstractDomain");
  }

  static void send_C_to_A(C_t& inp_C, A_t& out_A);
  static void send_C_to_B(C_t& inp_C, B_t& out_B);
};

/******************************************************************************/
/* a user does not need to look at the code below                             */
/******************************************************************************/
template <typename P>
struct Iterate_F {
  /*
   * A wrapper for calling f that uses
   */

  using A_t = typename P::A::t;
  using C_t = typename P::C::t;
  inline static A_t inp_A; /* A placeholder for the input */
  Iterate_F() {}

  /* Make this struct callable */
  void operator()(C_t& inp_C, C_t& out_C){
    /* convert inp:C_T -> inp:A_t */
    /* enhancement: in the future we can make this function depends on PRNG */
    /* that is defined within this struct. So, we can change the function easily */
    P::send_C_to_A(inp_A, inp_C);
    P::f(inp_A, out_C);
  }
};


template <typename Problem>
struct Iterate_G : Problem{
  /* accessing Problem::g will generate an error*/
  // Problem::g; /* Original iteration function */
  // using Problem::send_C_to_B;
  using B_t = typename Problem::B::t;
  using C_t = typename Problem::C::t;

  // todo should not be local to each thread?
  inline static B_t inp_B{}; /* A placeholder for the input */

  Iterate_G() {}

  /* Make this struct callable */
  void operator()(C_t& inp_C, C_t& out_C){

    /* convert inp:C_T -> inp:A_t */
    /* enhancement: in the future we can make this function depends on PRNG */
    /* that is defined within this struct. So, we can change the function easily */
    Problem::send_C_to_B(inp_B, inp_C);
    Problem::g(inp_B, out_C);
  }
};


template<typename C_t>
inline void swap_pointers(C_t*& pt1,
                          C_t*& pt2){
    /// pt1 will point to what pt2 was pointing at, and vice versa.
    C_t* tmp_pt = pt1;
    pt1 = pt2;
    pt2 = tmp_pt;
}

/* todo pass the function number. e.g. first version of (f, g) second version */
template <typename Problem>
void iterate_once(typename Problem::C::t* inp_pt,
		  typename Problem::C::t* out_pt,
		  Iterate_F<Problem>& F,
		  Iterate_G<Problem>& G)
{
  /*
   * Do 1 iteration inp =(f/g)=> out, write the output in the address pointed
   * by out_pt.
   */
  int f_or_g = Problem::C::extract_1_bit(*inp_pt);
  if (f_or_g == 1){
    F(*inp_pt, *out_pt);
    ++n_f_called;
  }
  else { /* f_or_g == 0 */
    G(*inp_pt, *out_pt);
    ++n_g_called;
  }
}




template<typename Problem>
auto generate_dist_point(typename Problem::C::t& inp0, /* don't change the pointer */
			 typename Problem::C::t*& tmp_inp_pt,
			 typename Problem::C::t*& out_pt,
                         u8* const output_bytes, /* size Problem::C::length */
			 u64& chain_length, /* write chain lenght here */
			 const i64 difficulty, /* #bits are zero at the beginning */
                         Iterate_F<Problem>& F,
                         Iterate_G<Problem>& G)
                         -> bool
{
  /*
   * Given an input, iterate functions either F or G until a distinguished point
   * is found, save the distinguished point in out_pt and output_bytes
   * Then return `true`. If the iterations limit is passed, returns `false`.
   */



  /* copy the input to tmp, then never touch the inp again! */
  Problem::C::copy(inp0, *tmp_inp_pt);

  static const u64 mask = (1LL<<difficulty) - 1;
  using C = typename  Problem::C;
  chain_length = 0;

  /* F: Domain_C -> Domain_C, instead of C_t -> A_t -f-> C_t */

  bool found_distinguished = false;

  found_distinguished =
    (0 == (mask & C::extract_k_bits(*tmp_inp_pt, difficulty))  );

  /* 1 if the next function is f */

  /* The probability of not finding a distinguished point during this loop is*/
  /* 2^{-some big k} */
  for (i64 i = 0; i < 40*(1LL<<difficulty); ++i){
    iterate_once(tmp_inp_pt, out_pt, F, G);

    ++chain_length;

    /* we may get a dist point here */
    found_distinguished
      = (0 == (mask & C::extract_k_bits(*tmp_inp_pt, difficulty))  );

    /* unlikely with high values of difficulty */
    if (found_distinguished) [[unlikely]]{
      C::serialize(*out_pt, output_bytes);
      return true; /* exit the whole function */
    }

    /* swap inp and out */
    swap_pointers(tmp_inp_pt, out_pt);
  }
  return false; /* no distinguished point were found */
}






template<typename Problem>
void walk(typename Problem::C::t*& inp1_pt,
	  typename Problem::C::t*& tmp1_pt, /* inp1 calculation buffer */
	  const u64 inp1_chain_len,
          typename Problem::C::t*& inp2_pt,
	  typename Problem::C::t*& tmp2_pt, /* inp2 calculation buffer */
	  const u64 inp2_chain_len,
          Iterate_F<Problem>& F,
          Iterate_G<Problem>& G)
{
  /* Given two inputs that lead to the same distinguished point,
   * find the earliest collision in the sequence before the distinguished point
   * add a drawing to illustrate this.
   */


  using C_t = typename Problem::C::t;


  size_t const diff_len = std::max(inp1_chain_len, inp2_chain_len)
                        - std::min(inp1_chain_len, inp2_chain_len);

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
  
 
  if (inp1_chain_len > inp2_chain_len){
    for (size_t i = 0; i < diff_len; ++i){
      iterate_once(inp1_pt, tmp1_pt, F, G);
      swap_pointers(inp1_pt, tmp1_pt);
    }
   }

 
  if (inp1_chain_len < inp2_chain_len){
    for (size_t i = 0; i < diff_len; ++i){
      iterate_once(inp2_pt, tmp2_pt, F, G);
      swap_pointers(inp1_pt, tmp1_pt);
    }
 
  }
  /****************************************************************************/
  
  /* now both inputs have equal amount of steps to reach a distinguished point */
  /* both sequences needs exactly `len` steps to reach dist point */
  size_t len = std::min(inp1_chain_len, inp2_chain_len);

  /* Check if the two inputs are equal then we have a robin hood */


  if(Problem::C::is_equal( *inp1_pt, *inp1_pt ))
    return; /* Robinhood */

  
  for (size_t i = 0; i < len; ++i){
    /* walk them together and check each time if their output are equal */
    /* return as soon equality is found */
    iterate_once(inp1_pt, tmp1_pt, F, G);
    iterate_once(inp2_pt, tmp2_pt, F, G);
    
    if(Problem::C::is_equal( *tmp1_pt, *tmp2_pt ))
      return; /* They are equal */
  }
}




template <typename Problem >
bool treat_collision(typename Problem::C::t*& inp1_pt,
		     typename Problem::C::t*& tmp1_pt, /* inp1 calculation buffer */
		     const u64 inp1_chain_len,
		     typename Problem::C::t*& inp2_pt,
		     typename Problem::C::t*& tmp2_pt, /* inp2 calculation buffer */
		     const u64 inp2_chain_len,
		     Iterate_F<Problem>& F,
		     Iterate_G<Problem>& G,
                     std::vector< std::pair<typename Problem::C::t,
		                            typename Problem::C::t> >& container)
{
  using A_t = typename Problem::A::t;
  using B_t = typename Problem::B::t;
  using C_t = typename Problem::C::t;


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
  /* walk inp1 and inp2 just before `x` */
  /* i.e. iterate_once(inp1) = iterate_once(inp2) */
  walk<Problem>(inp1_pt,
		tmp1_pt,
		inp1_chain_len,
		inp2_pt,
		tmp2_pt, /* inp2 calculation buffer */
		inp2_chain_len,
		F,
		G);
					   

  std::cout << "before pushback\n";


  /* assume copying */
  std::pair<C_t, C_t> p{*inp1_pt, *inp2_pt};
  
  /* todo here we should add more tests */
  container.push_back(std::move(p)); 
  return true;
}



/* ----------------------------- Naive method ------------------------------ */
template <typename T1, typename T2>
auto count_unique_2nd_elm(std::vector<std::pair<T1, T2>> arr) -> size_t
{
  if (arr.size() == 0) return 0;

  size_t n_unique_elm = 1;
  T2 tmp{arr[0].second};

  for (auto& elm_pair: arr ){
    if (tmp != elm_pair.second)
      ++n_unique_elm;
    tmp = elm_pair.second;
  }

  return n_unique_elm;
}

template <typename Problem, typename A_t, typename C_t,
          int length,               /* output length in bytes */
          size_t A_n_elements,
	  auto f, /* f: A_t -> C_t */
          auto next_A,
	  auto ith_elm> /* next_A(A_t& inp), edit inp to the next input  */

auto inp_out_ordered() /* return a list of all f inputs outputs */
  -> std::vector<std::pair<A_t, std::array<u8, length> > >

{
  /* this vector will hold all (inp, out) pairs of F */


  /* number of bytes needed to code an element of C */
  using C_serial = std::array<u8, length>;

  auto begin = wtime();
  std::vector< std::pair<A_t, C_serial> > inp_out(Problem::C::n_elements);


  auto end = wtime();
  auto elapsed_sec = end - begin;
  std::cout << "Initializing the first vector took " <<  elapsed_sec <<"sec\n";



  std::cout << "Starting to fill the A list "
            << A_n_elements
            << " elements\n";

  omp_set_num_threads(omp_get_max_threads());
  std::cout << "We are going to use " << omp_get_max_threads() << " threads\n";


  const int nthds = omp_get_max_threads();
  begin = wtime();
  #pragma omp parallel for schedule(static)
  for (int thd = 0; thd < nthds; ++thd){
    size_t offset = thd * (A_n_elements/nthds);
    size_t end_idx = (thd+1) * (A_n_elements/nthds);
    if (thd == nthds - 1) end_idx = A_n_elements;

    A_t inp{}; /* private variable for a thread */
    C_t out{}; /* private variable for a thread */

    ith_elm(inp, offset); /* set inp_A to the ith element */

    for (size_t i = offset; i < end_idx; ++i){

      Problem::f(inp, out);

      inp_out[i].first = inp;
      Problem::C::serialize(out, inp_out[i].second);


      /* get ready for the next round */
      next_A(inp);


    }
  }

  end = wtime();
  elapsed_sec  = end - begin;

  std::cout << std::fixed << "Done with the first list in time = " << elapsed_sec
            <<" sec\n";


  begin = wtime();
  /* sort the input output according to the second element (the output ) */
  std::sort(std::execution::par,
	    inp_out.begin(),
	    inp_out.end(),
	    [](std::pair<A_t, C_serial> io1, std::pair<A_t, C_serial> io2)
	    {return io1.second < io2.second; });

  end = wtime();
  elapsed_sec  = end - begin;
  std::cout << std::fixed << "Sorting took " << elapsed_sec << "sec\n";

  begin = wtime();
  size_t n_unique = count_unique_2nd_elm(inp_out);
  elapsed_sec  = wtime() - begin;
  std::cout << "counting #unique elements took " << elapsed_sec << "sec. it has "
	    << n_unique << " elements\n";



  /* Hopefully, NRVO will save the unwanted copying */
  return inp_out;

}









/* get all collisions by naive algorithm */
template <typename Problem>
auto all_collisions_by_list()
  -> std::vector< std::tuple<typename Problem::A::t,
			     typename Problem::B::t,
			     typename Problem::C::t>>

{
  /* this vector will hold all (inp, out) pairs of F */
  using A_t = typename Problem::A::t;
  using B_t = typename Problem::B::t;
  using C_t = typename Problem::C::t;


  /* number of bytes needed to code an element of C */
  const int nbytes = Problem::C::length;
  using C_serial = std::array<u8, nbytes>;

  /* collision container */

  std::vector< std::pair<A_t, C_serial> >
    inp_out_f_A_C = inp_out_ordered<Problem,
				    A_t,
				    C_t,
				    Problem::A::length,
				    Problem::A::n_elements,
				    Problem::f,
				    Problem::A::next,
				    Problem::A::ith_elm>();



  std::vector< std::pair<B_t, C_serial> >
    inp_out_g_B_C = inp_out_ordered<Problem,
				    B_t,
				    C_t,
				    Problem::B::length,
				    Problem::B::n_elements,
				    Problem::g,
				    Problem::B::next,
				    Problem::B::ith_elm>();


  std::cout << "Done createing the two lists!\n";
  /* now let's test all inputs of g and register those gets a collision */

  std::cout << "let's see the first 10 output sorted\n";
  for (int i = 0; i < 10; ++i){
    print_array(inp_out_f_A_C[i].second);
    print_array(inp_out_g_B_C[i].second);
    std::cout << "----------\n";
  }



  auto begin = wtime();
  auto end = wtime();
  auto elapsed_sec = end - begin;




  std::vector< std::tuple<A_t, B_t, C_t>  > all_collisions_vec{};
  C_t val{};
  size_t A_n_elements = inp_out_f_A_C.size();
  size_t B_n_elements = inp_out_g_B_C.size();

  begin = wtime();

  size_t idx_A = 0;
  size_t idx_B = 0;


  C_serial arr{};
  std::cout << "length of c output = " << arr.max_size() << ", or size = " << arr.size()
	    << "nbytes = " << nbytes << "\n" ;

  auto start = wtime();
  while ((idx_A < A_n_elements) and (idx_B < B_n_elements)) {
    /* 1st case: we have a collision  */
    if( inp_out_f_A_C[idx_A].second ==  inp_out_g_B_C[idx_B].second ) {
      std::cout << "Found collision at idx_A=" << idx_A << ", idx_B=" << idx_B
		<< " Do they collide? "<< (inp_out_f_A_C[idx_A].second ==  inp_out_g_B_C[idx_B].second)
		<<"\n";
      print_array(inp_out_f_A_C[idx_A].second);
      print_array(inp_out_g_B_C[idx_B].second);
      std::cout << "============================\n";

      Problem::C::unserialize(inp_out_f_A_C[idx_A].second, val);
      std::tuple col{inp_out_f_A_C[idx_A].first, inp_out_g_B_C[idx_B].first,  val };
      all_collisions_vec.push_back(col);

      ++idx_A;
      ++idx_B;
    }
    std::cout << "idx_A = " << idx_A << ", idx_B = " << idx_B << "\n";
    /* when */
    /* the list of inputs is sorted in increasing order */
    if( inp_out_f_A_C[idx_A].second <  inp_out_g_B_C[idx_A].second ){
      ++idx_A;
    } if( inp_out_f_A_C[idx_A].second >  inp_out_g_B_C[idx_A].second ){
      ++idx_B;
    }
  }

  elapsed_sec = wtime() - start;
  std::cout << std::fixed  << "it took " << elapsed_sec << " sec to find all collisions\n";
  std::cout << std::fixed <<"total " << all_collisions_vec.size() << " collisions\n";





  return all_collisions_vec;
}

/* --------------------------- end of Naive method ---------------------------- */



/* todo start from here */
template<typename Problem>
auto collision() -> std::pair<typename Problem::C::t, typename Problem::C::t>
{
  using A_t = typename Problem::A::t;
  using B_t = typename Problem::B::t;
  using C_t = typename Problem::C::t;
  /* Sanity Test: */
  std::cout << "unserial(serial(.)) =?= id(.) : "
	    << is_serialize_inverse_of_unserialize<Problem>()
	    << "\n";

  
  /************************* EXPERIMENT BY NAIVE METHOD ***********************/

  // std::cout << "Starting with naive method ...\n";
  // auto begin = wtime();
  // auto list_of_collisions = all_collisions_by_list<Problem>();
  // auto end = wtime();
  // auto elapsed = end - begin;

  // std::cout << std::fixed
  // 	    << "The naive method tells us that there are "
  //           << list_of_collisions.size()
  //           << " collisions. Took: "
  //           << elapsed
  //           << "s\n";

  /* end of EXPERIMENT BY NAIVE METHOD */
  /****************************************************************************/




  // --------------------------------- INIT -----------------------------------/
  // DICT
  /* todo get */

  // size_t n_slots = 1LL<<30; 

  size_t n_bytes = 1LL<<34; /* */
  Dict<u64, C_t> dict{n_bytes}; /* create a dictionary */



  // -----------------------------------------------------------------------------/
  // VARIABLES FOR GENERATING RANDOM DISTINGUISHED POINTS
  int difficulty = 2; // difficulty;
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
  /* Recall: Problem::C::length = #needed bytes to encode an element of C_t */
  u8 out0_bytes[Problem::C::length];


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


  /* Iteration Functions */
  Iterate_F<Problem> F{};
  Iterate_G<Problem> G{};
  /*********************************************************
   * surprise we're going actually use
   * void iterate_once(typename Problem::C::t* inp_pt,
   *		  typename Problem::C::t* out_pt,
   *		  Iterate_F<Problem>& F,
   *		  Iterate_G<Problem>& G)
   *
   *********************************************************/


  
  /*------------------- Generate Distinguished Points ------------------------*/
  while (n_collisions < n_needed_collisions){
    /* fill the input with a fresh random value. */
    // Problem::C::randomize(*inp_pt); 
    /* or use `/dev/urandom` */
    inp0 = read_urandom<C_t>();




    generate_dist_point<Problem>(inp0,
				 tmp0_pt,
				 out0_pt,
				 out0_bytes,
				 chain_length0,
				 difficulty,
				 F,
				 G);

    ++n_dist_points;
    //std::cout << "af inp = " << (*pt_inp_C)[0] << ", " << (*pt_inp_C)[1] << "\n";
    //std::cout << "bf out = " << (*pt_out_C)[0] << ", " << (*pt_out_C)[1] << "\n";
    //return std::pair(pt_inp_C, pt_inp_C);
    /* send the result to dictionary, check if it has a collision  */
    
    out0_digest = Problem::C::
      extract_k_bits(*out0_pt, Problem::C::length);
    
    is_collision_found = dict.pop_insert(out0_digest, /* key */
					 inp0, /* value  */
					 *inp1_pt, /* save popped element here */
					 chain_length1 /* of popped input */);
      
    /* todo work from here */
    
    if (is_collision_found) [[unlikely]]{
      ++n_collisions; 
      bool is_false_collision = false;

      /* respect the rule that inp0 doesn't have pointers dancing around it */
      Problem::C::copy(inp0, *tmp0_pt); /* (*tmp0_pt) holds the input value  */

      treat_collision<Problem>(tmp0_pt, /* inp0 */
			       out0_pt, /* inp0 scratch buffer  */
			       chain_length0,
			       inp1_pt,
			       out1_pt,
			       chain_length1,
			       F,
			       G,
			       collisions_container);


      std::cout << "found a collision " << n_collisions
		<<  " out of " << n_needed_collisions << "\n";

      std::cout << "false collisions = " << false_collisions
		<< " real collisions = " << real_collisions
		<< " robinhoods = " << n_robinhood
		<< "\n";

      std::cout << "#dist points = " << n_dist_points
		<<  " = 2^" << std::log2(n_dist_points) << "\n";

    }



  }


  return std::pair<C_t, C_t>(*tmp0_pt, *inp1_pt); // todo wrong values
}
#endif

// to use parallel sort sort
// install tbb lib
// sudo apt install libtbb-dev
// compiling
// g++ -flto -O3 -std=c++17  -fopenmp demos/speck32_demo.cpp -o speck32_demo -ltbb

// thd0 has out = 9ace,f6f8 and inp_out[i] = ce9a,f8f6 and inp = 1,0 and inp_out[i].first = 1,0
//AFTER unserialize thd0           has out = 9ace,f6f8 and inp_out[i] = ce9a,f8f6 and inp = 1,0 and inp_out[i].first = 1,0
