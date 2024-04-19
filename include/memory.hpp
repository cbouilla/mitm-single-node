/* get the available memory on a system. Currently only supports linux */
#ifndef MEMORY_INO
#define MEMORY_INO
#include <cstddef>
#include <linux/sysinfo.h>
#include <sys/sysinfo.h>
#include <cstdint>
#include <iostream>

inline /* get around  ODR (One Definition Rule) */
std::size_t get_available_memory()
{
  struct sysinfo info;
  sysinfo(&info);
  return info.freeram;
}


/* If nbytes exceeds the available memory, reduce it to available memeory. */
inline void adjust_to_available_memory(std::size_t &nbytes)
{
  if (nbytes > get_available_memory()){
    std::cout << "Adjusted memory from " << nbytes << " bytes to "
	      << get_available_memory() << " bytes\n";
    nbytes = get_available_memory();
  }
}

/* If nbytes was not specifed use all the available memory */
inline void set_to_available_memory(size_t& nbytes_memory)
{
  nbytes_memory = get_available_memory();
  printf("You have not specified how many bytes to be used for the dictionary!\n"
	 "Going to use %lu bytes\n" , nbytes_memory);

}

#endif 
