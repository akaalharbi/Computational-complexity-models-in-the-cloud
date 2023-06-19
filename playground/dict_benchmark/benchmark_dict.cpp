#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <bits/types/struct_timeval.h>
#include <stdint.h>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>
#include <unordered_map>
#include <cmath>

double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, 0);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}




int main ()
{
  // getrandom(void *buffer, size_t length, unsigned int flags)
  std::mt19937_64 rd;

  size_t  n = (1LL<<25);
  std::vector<uint64_t> inp(n);
  std::vector<uint64_t> qur(n);

  std::cout << "rd=" << rd() << "\n";
  double timing = wtime();
  
  for (size_t i = 0; i < inp.size(); ++i){
    inp[i] = rd();
    qur[i] = rd();
  }
  timing = wtime() - timing;
  std::cout << "It took "<< timing << "sec to fill random arrays\n";


  { /* default unordered_map */

    std::unordered_map<uint64_t, uint8_t> dict;

    timing = wtime();
    for (auto i = inp.begin(); i != inp.end(); ++i) {
      dict.emplace(std::make_pair(*i, 1));
    }
    timing = wtime() - timing;
    std::cout << "insert: std::unordered_map took " <<  timing << "sec \n";
    std::cout << "i.e. 2^" << std::log2(n/timing) << "elm/sec \n";

    uint64_t sum = 0;
    timing = wtime();
    for (auto i = qur.begin(); i != qur.end(); ++i) {
      auto found = dict.find(*i);
      if (found != dict.end()) 
	sum += found->second;
    }
    timing = wtime() - timing;
    std::cout << "lookup: std::unordered_map took " <<  timing << "sec \n";
    std::cout << "i.e. 2^" << std::log2(n/timing) << "elm/sec \n";
    std::cout << "sum = " << sum << "\n";
  }


  { /* https://probablydance.com/2017/02/26/i-wrote-the-fastest-hashtable/ */
    
  }

  { /* our dict :) */
    
  }
}

