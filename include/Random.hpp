#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

class Random
{
private:
  std::mt19937_64 re;
  std::uniform_real_distribution<double> urng;
public:
  Random(const long unsigned int seed)
  : re{seed}, urng{0.,1.} {}
  ~Random() {}
  inline double operator()() {return urng(re);}
};

#endif