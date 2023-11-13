/// Primitive adaptation of speck blockcipher using 16-bit words.
/// speck32* : {0, 1}^32 x {0, 1}^32 -> {0, 1}^32, the official documentation
/// specifies different length for k, we are using shorter key for simplicity
/// since security is not a concern.

#include <stdint.h>
#include <iostream>
#include <array>
#include "../mitm_sequential.hpp"
#include "speck32.hpp"



/******************************************************************************/
/* Setting up the problem                                                     */
/******************************************************************************/
using Speck32_t = typename std::array<uint16_t , 2>;

class SPECK_DOMAIN : AbstractDomain< Speck32_t >{

public:
    using t = Speck32_t;
    static const int length = 32;
    static const size_t n_elements = (1LL<<32);

    static bool is_equal(const t& x, const t& y){
       return x == y;
    }

    static void randomize(t& x) {
        x[0] = rand();
        x[1] = rand();
    }

    static void serialize(const t& x, uint8_t * out){
        out[0] = x[0];
        out[0] = x[0]>>8;
        out[1] = x[1];
        out[1] = x[1]>>8;
    }
    static void unserialize(t& out, const uint8_t* in){
        out[0] = in[0] | ((uint16_t ) in[1])<<8;
        out[1] = in[2] | ((uint16_t ) in[3])<<8;
    }

    inline static auto extract_1_bit(const t& inp) -> int {
        return inp[0]&1;
    }

    inline static auto extract_k_bits(const t& inp, int k) -> uint64_t {
        /* k = 16j + r, we would like to get the values of r and j  */
        k = k+1; /* Read bits after the first bit */
        uint16_t nbits_first_word = k&(16 - 1); /* read it mod 16 */
        /* first remove r and 16 at once by division, then make sure number < 16 */
        uint16_t nbits_second_word = (k>>4)&(16 - 1);

        /* What bits should we consider */
        uint16_t mask1 = (1<<nbits_first_word) - 1;
        uint16_t mask2 = (1<<nbits_second_word) - 1;

        /* maximally extrat 32 bits */
        return (inp[0]>>1)&mask1 | ((uint64_t) inp[1]&mask2)<<16;
    }
};

class Problem : AbstractProblem<SPECK_DOMAIN, SPECK_DOMAIN, SPECK_DOMAIN >{
public:

    /* Having the exact lines in the parent class doesn't help :( */
    using A = SPECK_DOMAIN;
    using A_t = typename A::t;

    using B = SPECK_DOMAIN;
    using B_t = typename B::t;

    using C = SPECK_DOMAIN;
    using C_t = typename C::t;

    static void f(const A_t &x, C_t &y){
        static std::array<uint16_t, 2> inpt_text{0, 0};
        encrypt(inpt_text, y, x);
    }

    static void g(const B_t &x, C_t &y){
        static std::array<uint16_t, 2> inp_ciphertext{0, 0};
        decrypt(inp_ciphertext, y, x);
    }

    static void send_C_to_A(A_t& out_A, C_t& inp_C){ out_A = inp_C; }
    static void send_C_to_B(B_t& out_B, C_t& inp_C){ out_B = inp_C; }
};


int main(int argc, char* argv[]){
    collision<Problem>();

}

