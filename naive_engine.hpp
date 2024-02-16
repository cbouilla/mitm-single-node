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

/* Generate all collisions/claws using naive algorithm. */



/* sorted search */ 

/* return Image(f) where f: X -> Y. */
template<typename X_t, typename Y_t>
auto images(void (*f)(X_t const& x, Y_t& y), /* y <- f(x) */
	    void (*ith_element)(X_t&, size_t),/* x = i */
	    void (*copy)(Y_t& y1, Y_t& y2), /* y2 := y1 */
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
      copy(y, outputs[ i + offset ]);
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

/* return -1 if a < b, 0 if a == b, 1 if a > b */
Ootemplate <typename Domain>
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

 
/* return {y: there is u, v satisfies y = f(u) = g(v) } sorted by the image where f: X -> Y. */
template <typename X_t, typename Y_t>
auto collisions(void (*copy)(const Y_t &y1, Y_t &y2),     /* y2 := y1 */
		int (*cmp)  (const Y_t& y1, const Y_t& y2),
		const std::vector<Y_t>& f_images,
		const std::vector<Y_t>& g_images		
		) -> std::vector<Y_t>
{
  std::vector<Y_t> collisions{};
  size_t n_elements_f = f_images.size();
  size_t n_found_collisions = 0;
  
   /* Parallelization could start from here */
  size_t nthreads = 1;
  for (size_t thd=0; thd < nthreads; ++thd){
    std::vector<Y_t> collisions_thd{};
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
				    [cmp](Y_t& a, Y_t& b){return (cmp(a, b) == -1);}); // cmp_with_start_f);

    auto end_g = std::upper_bound(g_images.begin(),
				  g_images.end(),
				  [cmp](Y_t& a, Y_t& b){return (cmp(a, b) < 0);});

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
      f_images = std::vector<Y_t>();
      g_images = std::vector<Y_t>();
      collisions.reserve(n_found_collisions);
    } /* thd0 is always responsible for the grand list of collisions */
    #pragma omp barrier /* all threads should've finished dealing with images_f, and images_g */
    

      
    /* merge vectors  */
    
   #pragma omp critical
    {
      collisions.insert(collisions.begin(), collisions_thd.begin(), collisions_thd.end());
    }
  } /* end of threads work */
  /* sort collisions in parallel */
  std::sort(std::execution::par,
	    collisions.begin(),
	    collisions.end(),
	    [cmp](Y_t const& a, Y_t const& b){return (cmp(a, b) < 0);}
	    );

  return collisions;
}


/*
 * Given a problem that follows AbstractDomain, and AbstractCollisionProblem.
 * Try to find two inputs x and y s.t. f(x) = f(y), where:
 * f: A -> C using a naive algortihm
 * Note: It's guaranteed in this method to find all collisions.
 */
template <typename Problem>
void naive_collisoin_search(Problem& Pb)
{
  /* aïe aïe aïe, just realized that it won't work looking at the inputs! */
}


/*
 * Given a problem that follows AbstractDomain, and AbstractClawProblem.
 * Try to find two inputs x_A and x_B s.t. f(x_A) = g(x_B), i.e.
 * a claw between f and g.
 * Note: It's guaranteed in this method to find all claws.
 */
template <typename Problem> void claw_search(Problem &Pb)
{
  /* get all the images of f */
  /* get all the images of g */
  /* compare them */
  /* save the results in disk! */
}




#endif




