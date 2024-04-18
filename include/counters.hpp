#ifndef MITM_COUNTERS
#define MITM_COUNTERS
#include "timing.hpp"
#include "util_file_system.hpp"
#include <numeric>
#include <vector>
#include <fstream>

namespace mitm {
/* Non-essential counters but helpful to have, e.g. n_collisions/sec */
struct Counters {
  size_t const interval = (1LL<<15); // for printing only
  size_t n_updates = 0; // #times dict were flushed
  size_t n_collisions = 0;
  
  // each entry contains number of distinguished point between
  size_t n_dist_points_previous = 0;
  std::vector<size_t> n_distinguished_points = {0};

  double start_time;
  double end_time;
  double dist_previous_time;
  double update_previous_time;
  double elapsed = 0;

  Counters()
      : start_time(wtime()), dist_previous_time(wtime()),
        update_previous_time(wtime())
        {}

  Counters(double interval)
    : interval(interval),
      dist_previous_time(wtime()),
      update_previous_time(wtime())
  {}


  void increment_n_distinguished_points()
  {
    ++n_distinguished_points[n_updates];
    size_t n = n_distinguished_points[n_updates] - n_dist_points_previous;

    if (n_distinguished_points[n_updates] % interval == 0){

      elapsed = wtime() - dist_previous_time;
      printf("\r%lu=2^%0.2f iter took:"
	     " %0.2f sec, i.e. %0.2f ≈ 2^%0.2f iter/sec",
	     n, std::log2(n),
	     elapsed, n/elapsed, std::log2(n/elapsed) );

      fflush(stdout);
      n_dist_points_previous = n_distinguished_points[n_updates];
      dist_previous_time = wtime();
    }
  }

  
  void increment_n_updates()
  {
    elapsed = wtime() - update_previous_time;

    /* new entry to */
    printf("\nUpdating the iteration function.\n"
	   "the previous iteration:\n"
	   " - lasted %0.2f sec\n"
	   " - generated %lu ≈ 2^%0.2f distinguished points\n",
	   elapsed,
	   n_distinguished_points[n_updates],
	   std::log2(n_distinguished_points[n_updates]));
    
    ++n_updates;
    n_distinguished_points.emplace_back(0);
    update_previous_time = wtime();
  }


  void increment_collisions(size_t n = 1) {n_collisions += n;}

  
  void save_summary_stats(std::string problem_type,
			  size_t A_size,
			  size_t B_size,
			  size_t C_size,
			  int difficulty)
  {
    end_time = wtime();
    double total_time = end_time - start_time;
    printf("----------------------------------------\n"
	   "Took %0.2f sec to find the golden inputs.\n"
	   "Saving the counters ...\n",
	   total_time);
    /* open folder (create it if it doesn't exist )*/
    std::ofstream summary;
    std::string d_name = "data/";
    std::string f_name = d_name + problem_type + "_summary.csv";
    create_folder_if_not_exist(d_name);
    int file_status = create_file_if_not_exist(f_name);

      
    /* open the summary file */
    summary.open(f_name, std::ios::app);


    std::string column_names = "";
    /* Depending on the problem, we have different column names */
    if (problem_type == "claw")
      column_names = "C_size,A_size,B_size,difficulty,#distinguished_points,log2(#distinguished_points),#collisions,log2(#collisions),#updates,time(sec)\n";
    if (problem_type == "collision")
      column_names = "C_size,A_size,difficulty,#distinguished_points,log2(#distinguished_points),#collisions,log2(#collisions),#updates,time(sec)\n";

    /* Write column names to the file only if the file did not exist before */
    if (file_status == 2)
      summary << column_names;


    /* compute some stats */
    size_t total_distinguished_points
      = std::accumulate(n_distinguished_points.begin(),
			n_distinguished_points.end(),
			static_cast<size_t>(0)); // starting value.

    double log2_n_distinguished_points = std::log2(total_distinguished_points);
    double log2_n_collisions = std::log2(n_collisions);

    /* write all those stats to the summary file */
    std::string B_data = "";
    if (B_size != 0) /* B_size passed as 0 when this function is called by a
		      * collision problem */
      B_data = std::to_string(B_size) +  ", ";
    
    summary << std::to_string(C_size) << ", "
	    << std::to_string(A_size) << ", "
	    << B_data
	    << std::to_string(difficulty) << ", "
	    << std::to_string(total_distinguished_points) << ", "
      	    << std::to_string(log2_n_distinguished_points) << ", "
      	    << std::to_string(n_collisions) << ", "
	    << std::to_string(log2_n_collisions) << ", "
      	    << std::to_string(n_updates) << ", "
	    << total_time
	    << "\n";


    summary.close();

    std::cout << "Successfully saved stats in " << f_name << "\n"
	      << "Format:\n"
	      << column_names
	      << "\n" /* This should be the end */
	      << std::to_string(C_size) << ", "
	      << std::to_string(A_size) << ", "
	      << std::to_string(B_size) << ", "
	      << std::to_string(difficulty) << ", "
	      << std::to_string(total_distinguished_points) << ", "
	      << std::to_string(log2_n_distinguished_points) << ", "
	      << std::to_string(n_collisions) << ", "
	      << std::to_string(n_updates) << ", "
	      << total_time
	      << "\n";
  }
};  
}
#endif
