#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>

class Random
{
private:
  std::mt19937_64 re;
  std::uniform_real_distribution<double> urng;
  std::uniform_int_distribution<int> iurng;
  int calls;
public:
  Random(const long unsigned int seed)
  : re{seed}, urng{0.,1.}, iurng{1,5}, calls{0} {}
  ~Random() {}
  inline double operator()() {
    calls +=1;
    return urng(re);
  }
  inline int randint(){
    return iurng(re);
  }

  int GetCalls() const {return calls;}
};

#endif