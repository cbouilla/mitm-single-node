#ifndef NAIVE_ENGINE
#define NAIVE_ENGINE
#include <cstddef>
#include <vector>
#include <tuple>
#include "AbstractDomain.hpp"

/* Generate all collisions/claws using naive algorithm. */



/* sorted search */ 

/* return Image(f) where f: X -> Y. */
template<typename X_t, typename Y_t>
auto images(void (*f)(X_t const& x, Y_t& y), /* y <- f(x) */
	    void (*ith_element)(X_t&, size_t),/* x = i */
	    void (*copy)(Y_t& y1, Y_t& y2), /* y2 := y1 */
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
    size_t end = (n_elements/nthreads) + (n_elements%n_elements)*(thd == (nthreads -1));

    /* A thread main work */
    for (size_t i = 0; i < end; ++i){
      f(x, y);
      /* todo make it by index */
      copy(y, outputs[ i + offset ]);
      ith_element(x, i + global_offset );
    }
  }

  return outputs;
}


/* return {y: there is u, v satisfies y = f(u) = g(v) } sorted by the image where f: X -> Y. */
template <typename X_t, typename Y_t>
auto collisions(void (*g)(X_t const &x, Y_t &y),    /* y <- f(x) */
		void (*ith_element)(X_t &, size_t), /* x = i */
		void (*copy)(Y_t &y1, Y_t &y2),     /* y2 := y1 */
		size_t start,                       /* index of first x to treat */
		size_t end /* index of last x to treat. it's can be generated*/,
		std::vector<Y_t>& f_images
		) -> std::vector<Y_t>
{
  size_t n_elements = end - start; /* to be processed in this function */
  std::vector<Y_t> collisions{};
  


  /* Parallelization could start from here */
  size_t nthreads = 1;
  
  for (size_t thd = 0; thd < nthreads; ++thd){
    Y_t y{};
    X_t x{};
    std::vector<Y_t> collisions_prv{};
    const size_t offset = thd*(n_elements/nthreads);
    /* get the ith input */
    ith_element(x, 0 + offset);
    /* last term accounts for r in  n_elements = nthreads*k + r, in the last thread */
    size_t end = (n_elements/nthreads) + (n_elements%n_elements)*(thd == (nthreads -1));

    /* A thread main work */
    for (size_t i = 0; i < end; ++i){
      g(x, y);
      /* todo add it to a local container  */
      // if found collisions using sorted search
      
      // merge collisions
      ith_element(x, i + offset );
    }
    /* merge local container with the global one */
    // 
  }

  return collisions;
}



#endif
