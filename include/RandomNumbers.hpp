#include <random>
#include <optional>

#ifndef RANDOMNUMBERS_H
#define RANDOMNUMBERS_H

using namespace std;

class RandomNumbers {
private:
  double range_lo;
  double range_hi;

  random_device *rd;
  mt19937 *generator;
  uniform_real_distribution<double> *uniform_dist;
  uint32_t seed;          // Stores the seed value
  bool is_seed_set;       // True if seed is set

public:
  RandomNumbers() {
    range_lo = 0.0;
    range_hi = 1.0;

    rd = nullptr;
    generator = nullptr;
    uniform_dist = nullptr;
    seed = 0;
    is_seed_set = false;
  }
  RandomNumbers(double _range_lo, double _range_hi) {
    range_lo = _range_lo;
    range_hi = _range_hi;

    rd = new random_device();
    generator = new mt19937((*rd)());
    uniform_dist = new uniform_real_distribution<double>(range_lo, range_hi);
    seed = 0;
    is_seed_set = false;
  }
  RandomNumbers(double _range_lo, double _range_hi, uint32_t _seed) {
    range_lo = _range_lo;
    range_hi = _range_hi;

    rd = nullptr;
    generator = new mt19937(_seed);
    uniform_dist = new uniform_real_distribution<double>(range_lo, range_hi);
    seed = _seed;
    is_seed_set = true;
  }
  ~RandomNumbers() {
    delete rd;
    delete generator;
    delete uniform_dist;
  }

  double Uniform() { return (*uniform_dist)(*generator); }

  uint32_t GetSeed() const { return seed; }
  bool IsSeedSet() const { return is_seed_set; }
};

#endif // RANDOMNUMBERS_H
