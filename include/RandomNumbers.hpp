#include <random>

#ifndef RANDOMNUMBERS_H
#define RANDOMNUMBERS_H

using namespace std;

class RandomNumbers {
private:
  double range_lo;
  double range_hi;

  random_device *rd;
  // mt19937 generator(rd());
  mt19937 *generator;
  // uniform_real_distribution<double> uniform_dist(range_lo, range_hi);
  uniform_real_distribution<double> *uniform_dist;

public:
  RandomNumbers() {
    /*
        Default Constructor
    */
    range_lo = 0.0;
    range_hi = 1.0;

    rd = nullptr;
    generator = nullptr;
    uniform_dist = nullptr;
  }
  RandomNumbers(double _range_lo, double _range_hi) {
    /*
        Constructor without seed, by default seeds from random device (standard
       way)
    */
    range_lo = _range_lo;
    range_hi = _range_hi;

    rd = new random_device();
    generator = new mt19937((*rd)());
    uniform_dist = new uniform_real_distribution<double>(range_lo, range_hi);
  }
  RandomNumbers(double _range_lo, double _range_hi, uint32_t _seed) {
    /*
        Constructor with seed, for reproducibility
    */
    range_lo = _range_lo;
    range_hi = _range_hi;

    rd = nullptr;
    generator = new mt19937(_seed);
    uniform_dist = new uniform_real_distribution<double>(range_lo, range_hi);
  }
  ~RandomNumbers() {
    delete rd;
    delete generator;
    delete uniform_dist;
  }

  double Uniform() { return (*uniform_dist)(*generator); }
};

#endif // RANDOMNUMBERS_H
