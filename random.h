#ifndef RANDOM_H
#define RANDOM_H

#include <random>
#include <algorithm>
#include "./SFMT-src-1.5.1/SFMT.h"

using namespace std;

class Random{
public:
  sfmt_t sfmt;
  Random(){
    sfmt_init_gen_rand(&sfmt, time(0));
  }
  unsigned int generateInt() {
    return sfmt_genrand_uint32(&sfmt);
  }
  double generateReal() {
    return sfmt_genrand_real2(&sfmt);
  }
};

#endif
