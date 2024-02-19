#ifndef NAIVE_ENGINE
#define NAIVE_ENGINE
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <vector>
#include <tuple>
#include <algorithm>
#include <execution>
#include <omp.h>
#include "AbstractDomain.hpp"


/* return -1 if a < b, 0 if a == b, 1 if a > b */
template <typename Domain>
int compare(const typename Domain::t& a,
	    const typename Domain::t& b)
{
  Domain dom{};
  static uint8_t a_serial[dom.length];
  static uint8_t b_serial[dom.length];

  dom.serialize(a, a_serial);
  dom.serialize(b, b_serial);

  return std::memcmp(a_serial, b_serial, dom.length);
}


/* given (a0, a1) and (b0, b1) return -1 if a0 < b0, 0 if a0 == b0, 1 if a0 > b0 */
template <typename Domain, typename Y_t, typename X_t>
int compare_pair(const typename std::pair<Y_t, X_t>& a,
		 const typename std::pair<Y_t, X_t>& b)
{/* todo take a1, b1 into considerations */
  Domain dom{};
  static uint8_t a_serial[dom.length];
  static uint8_t b_serial[dom.length];

  dom.serialize(a.first, a_serial);
  dom.serialize(b.first, b_serial);

  return std::memcmp(a_serial, b_serial, dom.length);
}


/* Generate all collisions/claws using naive algorithm. */



/* sorted search */ 

/* return Image(f) where f: X -> Y. */
template<typename X_t, typename Y_t>
auto images(void (*f)(X_t const& x, Y_t& y), /* y <- f(x) */
	    void (*ith_element)(X_t&, size_t),/* x = i */
	    bool (*cmp)(Y_t& y1, Y_t& y2),
	    size_t start, /* index of first x to treat */
	    size_t end /* index of last x to treat. it's can be generated*/
	    ) -> std::vector<Y_t>
{
  size_t n_elements = end - start; /* to be processed in this function */
  std::vector<Y_t> outputs{};
  outputs.resize(n_elements);


  /* Parallelization could start from here */
  size_t nthreads = 1;
  
  for (size_t thd = 0; thd < nthreads; ++thd){
    Y_t y{};
    X_t x{};

    /* thd j, will work on indices that starts with offset */
    const size_t offset =  thd*(n_elements/nthreads);
    const size_t global_offset = start + offset; /*  globally which x to treat */

    /* get the ith input */
    ith_element(x, 0 + global_offset);
    /* last term accounts for r in  n_elements = nthreads*k + r, in the last thread */
    size_t nelement_per_thd = (n_elements/nthreads) + (n_elements%n_elements)*(thd == (nthreads -1));

    /* A thread main work */
    for (size_t i = 0; i < nelement_per_thd; ++i){
      f(x, y);
      /* todo make it by index */
      outputs[ i + offset ] = y;
      ith_element(x, i + global_offset );
    }
  }

  /* sort outputs in parallel */
  std::sort(std::execution::par,
	    outputs.begin(),
	    outputs.end(),
	    cmp);
  

  return outputs;
}


/* returns (X, f(X))*/
template<typename X_t, typename Y_t>
auto domain_images(void (*f)(X_t const& x, Y_t& y), /* y <- f(x) */
		   void (*ith_element)(X_t&, size_t),/* x = i */
		   bool (*cmp)(Y_t& y1, Y_t& y2),
		   size_t start, /* index of first x to treat */
		   size_t end /* index of last x to treat. it's can be generated*/
		   ) -> std::vector< std::pair<X_t, Y_t> >
{
  size_t n_elements = end - start; /* to be processed in this function */
  std::vector< std::pair<X_t, Y_t> > outputs{};
  outputs.resize(n_elements);


  /* Parallelization could start from here */
  size_t nthreads = 1;
  
  for (size_t thd = 0; thd < nthreads; ++thd){
    Y_t y{};
    X_t x{};

    /* thd j, will work on indices that starts with offset */
    const size_t offset =  thd*(n_elements/nthreads);
    const size_t global_offset = start + offset; /*  globally which x to treat */

    /* get the ith input */
    ith_element(x, 0 + global_offset);
    /* last term accounts for r in  n_elements = nthreads*k + r, in the last thread */
    size_t nelement_per_thd = (n_elements/nthreads) + (n_elements%n_elements)*(thd == (nthreads -1));

    /* A thread main work */
    for (size_t i = 0; i < nelement_per_thd; ++i){
      f(x, y);
      /* todo make it by index */
      outputs[ i + offset ] = std::pair(y, x);
      ith_element(x, i + global_offset );
    }
  }

  /* sort outputs in parallel */
  std::sort(std::execution::par,
	    outputs.begin(),
	    outputs.end(),
	    cmp);
  

  return outputs;
}


 
/* return {y: there is u, v satisfies y = f(u) = g(v) } sorted by the image where f: X -> Y. */
//template <typename X_t, typename Y_t>
template <typename T>
auto extract_collisions(int (*cmp)  (const T& y1, const T& y2),
			const std::vector<T>& f_images,
			const std::vector<T>& g_images		
			) -> std::vector<T>
{
  std::vector<T> collisions{};
  size_t n_elements_f = f_images.size();
  size_t n_found_collisions = 0;
  
   /* Parallelization could start from here */
  size_t nthreads = 1;
  for (size_t thd=0; thd < nthreads; ++thd){
    std::vector<T> collisions_thd{};
    int result = 0;

    /* which indices of f to read: */
    size_t idx_start_f = thd*(n_elements_f/nthreads);
    size_t idx_end_f = thd*(n_elements_f/nthreads) + (n_elements_f%n_elements_f)*(thd == (nthreads -1));
    /* what indices of g to read */
    // constexpr ForwardIt lower_bound( ForwardIt first, ForwardIt last,
    //                                  const T& value, Compare comp );
    auto start_g = std::lower_bound(g_images.begin(),
				    g_images.end(),
				    f_images[idx_start_f],
				    [cmp](T& a, T& b){return (cmp(a, b) == -1);}); // cmp_with_start_f);

    auto end_g = std::upper_bound(g_images.begin(),
				  g_images.end(),
				  [cmp](T& a, T& b){return (cmp(a, b) < 0);});

    auto start_f = (f_images.begin() + idx_start_f);
    auto end_f = (f_images.begin() + idx_end_f);

    while ( (start_g != end_g) || (start_f != end_f) ) {
      /* if they are equal, store the collision */
      result = cmp(*start_f, *start_g);
      if (result == 0){
	collisions_thd.push_back(*start_f);
	++start_f;
	++start_g;
      }

      /* elif f < g, forward only f */
      if (result < 0)
	++start_f;

      /* elif g < f, forward only g */
      if (result > 0)
	++start_g;
    }


    #pragma omp atomic
    { n_found_collisions += collisions_thd.size(); }

    #pragma omp barrier /* all threads should've finished dealing with images_f, and images_g */
    if (thd == 0){
      /* free spaces used by images_g and images_f */
      f_images = std::vector<T>();
      g_images = std::vector<T>();
      collisions.reserve(n_found_collisions);
    } /* thd0 is always responsible for the grand list of collisions */
    #pragma omp barrier /* all threads should've finished dealing with images_f, and images_g */
    

      
    /* merge vectors  */
    
    #pragma omp critical
    {/* add thd found collisions to the shared list */
      collisions.insert(collisions.begin(),
			collisions_thd.begin(),
			collisions_thd.end());
    }
  } /* end of threads work */
  /* sort collisions in parallel */
  std::sort(std::execution::par,
	    collisions.begin(),
	    collisions.end(),
	    [cmp](T const& a, T const& b){return (cmp(a, b) < 0);}
	    );

  return collisions;
}


/*
 * Given a problem that follows AbstractDomain, and AbstractCollisionProblem.
 * Try to find two inputs x and y s.t. f(x) = f(y), where:
 * f: A -> C using a naive algortihm
 * Note: It's guaranteed in this method to find all collisions.
 * This function asks for an extra input "ith_element"
 * s.t. ith_element(x, i), x <- ith_input (according to some order, doesn't matter which )
 *      ith_element(x, i) =/= ith_element(x, j) if i =/= j.
 */
template <typename Problem>
auto naive_collisoin_search(Problem& Pb,
			    void (*ith_element)(typename Problem::A_t &, size_t),
			    size_t n_elements
			    ) -> std::vector<std::pair<typename Problem::C_t,
						       typename Problem::A_t>>
{

  using A_t = typename Problem::A_t;
  using C_t = typename Problem::C_t;
  using Dom_C = typename Problem::Dom_C;


  /* create all possible pairs (f(x), x) sorted by f(x) as std::vector<std::pair<>> */
  auto inps_outs = domain_images(Pb.f,
				 ith_element,
				 compare_pair<Dom_C, C_t, A_t>,
				 0, /* start, this argument for future use with MPI */
				 n_elements);

  /* return all pairs (f(x), x) s.t. there is an x' =/= x and f(x') == f(x), do this for all x */
  auto collisions = extract_collisions(compare_pair<Dom_C, C_t, A_t>, inps_outs, inps_outs);

  /* save collisions in disk */
  // todo
  return collisions;
}


/*
 * Given a problem that follows AbstractDomain, and AbstractClawProblem.
 * Try to find two inputs x_A and x_B s.t. f(x_A) = g(x_B), i.e.
 * a claw between f and g.
 * Note: It's guaranteed in this method to find all claws.
 */
template <typename Problem>
auto claw_search(Problem &Pb,
		 void (*ith_element)(typename Problem::A_t &, size_t),
		 size_t n_elements_A,
		 size_t n_elements_B
		 ) -> std::vector<typename Problem::C_t>
{

  using A_t = typename Problem::A_t;
  using B_t = typename Problem::B_t;
  using C_t = typename Problem::C_t;
  using Dom_C = typename Problem::Dom_C;


  /* get all f(x) where x in A, stored in std::vector<C_t>  */
  auto f_images = images(Pb.f,
			 ith_element,
			 compare<Dom_C>,
			 0, /* start with the 1st element in A */
			 n_elements_A); /* end with the last element in A */

  /* get all f(x) where x in A, stored in std::vector<C_t>  */
  auto g_images = images(Pb.g,
			 ith_element,
			 compare<Dom_C>,
			 0, /* start with the 1st element in A */
			 n_elements_B); /* end with the last element in A */

  
  


  /* todo iterating only over A inputs, may not return all claws!
   *  we should iterate over the largest domain
   */
  auto collisions = extract_collisions(compare<Dom_C>,
				       f_images,
				       g_images);

  /* todo save the results in disk! */


  return collisions;
  
}




#endif




