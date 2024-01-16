#ifndef MITM_PRNG
#define MITM_PRNG
#include <chrono>
#include <fstream>
#include <array>
#include <cstddef>
#include <cstdint>
#include<cstring>


namespace mitm {
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



/*
 * Generic interface for a PRNG. The sequence of pseudo-random numbers
 * depends on both seed and seq
 */
class PRNG {
  /* todo make it support length */
  /* source : https://gist.github.com/mortenpi/9745042 */
  union {
    uint64_t value;
    char cs[sizeof(uint64_t)];
  } u;
  
public:

  PRNG() { };
  
  uint64_t rand(){
    std::ifstream rfin("/dev/urandom");
    rfin.read(u.cs, sizeof(u.cs));
    rfin.close();

    return u.value;
  };
  
};


}

#endif
