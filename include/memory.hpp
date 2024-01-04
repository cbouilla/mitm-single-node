/* get the available memory on a system. Currently only supports linux */
#ifndef MEMORY_INO
#define MEMORY_INO
#include <linux/sysinfo.h>
#include <sys/sysinfo.h>
#include <cstdint>


static std::size_t get_available_memory()
{
  struct sysinfo info;
  sysinfo(&info);
  return info.freeram;
}

#endif 
