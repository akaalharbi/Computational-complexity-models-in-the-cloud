// g++ -O3 -std=c++20 benchmark_dict.cpp  -o benchmark_dict && ./benchmark_dict 
// g++ -march=native -flto -O3 -std=c++20 dict.c benchmark_dict.cpp  -o benchmark_dict && ./benchmark_dict
/* ---------------- RESULTS ----------------------- */
// raw
// ahmed@ahmed:~/phd_shared/playground/dict_benchmark$  g++ -march=native -flto
// -O3 -std=c++20 dict.c benchmark_dict.cpp  -o benchmark_dict &&
// ./benchmark_dict Memory usage: ram: 1716 kb vm: 6060 kb size of insert
//
// n = 2^28
// 1073741824 bytes
// ================================================
// after the array init
// Memory usage: ram: 2099604 kb vm: 2103220 kb
// ================================================
// It took 0.501639sec to fill random arrays
//
// --------------------------------------------
// - insert: our dict took 12.9237sec 
// i.e. 2^23.308elm/sec 
// - lookup: our dict took 1.51491sec 
// i.e. 2^26.4008elm/sec 
// sum = 1
// Memory usage: ram: 2624116 kb vm: 4200376 kb
// ================================================
// after exiting the dictionary
// Memory usage: ram: 2099844 kb vm: 2103220 kb
// ================================================

// ----------------------------------------------
// - insert: std::unordered_map took 48.1389sec
// i.e. 2^21.4109elm/sec
// - lookup: std::unordered_map took 9.73339sec
// i.e. 2^23.7171elm/sec
// size = 768614336404564650 bytes
// sum = 0
//  std::unorderd_map size is taken into account
// Memory usage: ram: 7903760 kb vm: 7907388 kb
//
// ================================================
// after exiting the dictionary
// Memory usage: ram: 6294148 kb vm: 6297524 kb
// ================================================
// ------------------------------------------------
//  // this person claims his hash is faster than google dense_hash
// - insert: ska::unordered_map took 28.1248sec 
// i.e. 2^22.1862elm/sec  // it reaches sometimes to 2^23 
// - lookup: ska::flat_hash_map took 6.27012sec 
// i.e. 2^24.3515elm/sec 
// - size = 32025597350190193 bytes
// sum = 0
//  ska::flat_hasp_map size is taken into account 
// Memory usage: ram: 12585612 kb vm: 12588984 kb
// ================================================
// after exiting the dictionary
// Memory usage: ram: 6294152 kb vm: 6297524 kb
// ================================================

//------------------------------------------------
// Trying to insert
// 1- std::unorderd_map
// n=2^30 fail
//
// 2- ska::flat_hash_map
//
// 3- our dict
// n = 2^30, u64
//- insert: our dict took 102.996sec
// i.e. 2^23.3136elm/sec
// - lookup: our dict took 66.5931sec
// i.e. 2^23.9427elm/sec
// n = 2^29, u64
// - insert: our dict took 44.4027sec 
// i.e. 2^23.5274elm/sec 
// - lookup: our dict took 33.4775sec 
// i.e. 2^23.9349elm/sec 
// n=28
// - insert: our dict took 16.9949sec 
// i.e. 2^23.913elm/sec 
// - lookup: our dict took 9.45432sec 
// i.e. 2^24.759elm/sec 

/* -------------------------------------------------- */


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
#include "flat_hash_map.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "dict.h"

double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, 0);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}


/// source: https://hpcf.umbc.edu/general-productivity/checking-memory-usage/
/*
 * Look for lines in the procfile contents like: 
 * VmRSS:         5560 kB
 * VmSize:         5560 kB
 *
 * Grab the number between the whitespace and the "kB"
 * If 1 is returned in the end, there was a serious problem 
 * (we could not find one of the memory usages)
 */
int get_memory_usage_kb(long* vmrss_kb, long* vmsize_kb)
{
  /* Get the the current process' status file from the proc filesystem */
  FILE* procfile = fopen("/proc/self/status", "r");

  long to_read = 8192;
  char buffer[to_read];
  int read = fread(buffer, sizeof(char), to_read, procfile);
  ++read; // dummy operation to avoid non used variable warning
  fclose(procfile);

  short found_vmrss = 0;
  short found_vmsize = 0;
  char* search_result;

  /* Look through proc status contents line by line */
  char delims[] = "\n";
  char* line = strtok(buffer, delims);

  while (line != NULL && (found_vmrss == 0 || found_vmsize == 0) )
    {
      search_result = strstr(line, "VmRSS:");
      if (search_result != NULL)
        {
	  sscanf(line, "%*s %ld", vmrss_kb);
	  found_vmrss = 1;
        }

      search_result = strstr(line, "VmSize:");
      if (search_result != NULL)
        {
	  sscanf(line, "%*s %ld", vmsize_kb);
	  found_vmsize = 1;
        }

      line = strtok(NULL, delims);
    }

  return (found_vmrss == 1 && found_vmsize == 1) ? 0 : 1;
}



void print_memory_usage(){
  
  long vmrss_kb, vmsize_kb;
  get_memory_usage_kb(&vmrss_kb, &vmsize_kb);
  printf("Memory usage: ram: %ld kb vm: %ld kb\n",
	  vmrss_kb, vmsize_kb);

}

int main ()
{
  // getrandom(void *buffer, size_t length, unsigned int flags)
  std::mt19937_64 mt;
  std::random_device rd;
  mt.seed(rd());
  
  
  print_memory_usage();
  size_t  n = (1LL<<27);
  std::vector<uint64_t> inp(n);
  std::vector<uint64_t> qur(n);

  std::cout << "size of insert " << inp.size()*sizeof(uint64_t) << " bytes\n";
  double timing = wtime();
  
  for (size_t i = 0; i < inp.size(); ++i){
    inp[i] = mt();
    qur[i] = mt();
  }
  timing = wtime() - timing;
  std::cout << "================================================\n" ;
  std::cout << "after the array init \n";
  print_memory_usage();
  std::cout << "================================================\n" ;

  std::cout << "It took "<< timing << "sec to fill random arrays\n";


  { /* our dictionary  */
    dict* d = dict_new(n);
    uint64_t sum = 0;
    timing = wtime();
    for (size_t i = 0; i<n; ++i) {
      dict_add_element_to(d, (u8*) &inp[i]);
    }
    timing = wtime() - timing;
    std::cout << "- insert: our dict took " <<  timing << "sec \n";
    std::cout << "i.e. 2^" << std::log2(n/timing) << "elm/sec \n";

    for (size_t i = 0; i<n; ++i) {
      d->values[i] = mt();
    }

    timing = wtime();
    for (size_t i = 0; i<n; ++i) {
      sum += dict_has_elm(d, (u8*) &qur[i]);
    }
    timing = wtime() - timing;
    std::cout << "- lookup: our dict took " <<  timing << "sec \n";
    std::cout << "i.e. 2^" << std::log2(n/timing) << "elm/sec \n";
    std::cout << "sum = " << sum << "\n";
    print_memory_usage();
    free(d->values);
    free(d);
    
   }

  std::cout << "================================================\n" ;
  std::cout << "after exiting the dictionary\n";
  print_memory_usage();
  std::cout << "================================================\n" ;
  

  { /* default unordered_map */

    std::unordered_map<uint64_t, uint8_t> dictionary;

    timing = wtime();
    for (auto i = inp.begin(); i != inp.end(); ++i) {
      dictionary.emplace(std::make_pair(*i, 1));
    }
    timing = wtime() - timing;
    std::cout << "- insert: std::unordered_map took " <<  timing << "sec \n";
    std::cout << "i.e. 2^" << std::log2(n/timing) << "elm/sec \n";

    uint64_t sum = 0;
    timing = wtime();
    for (auto i = qur.begin(); i != qur.end(); ++i) {
      auto found = dictionary.find(*i);
      if (found != dictionary.end()) 
	sum += found->first;
    }
    timing = wtime() - timing;
    std::cout << "- lookup: std::unordered_map took " <<  timing << "sec \n";
    std::cout << "i.e. 2^" << std::log2(n/timing) << "elm/sec \n";
    std::cout << "size = " << dictionary.max_size() << " bytes\n";
    std::cout << "sum = " << sum << "\n";
    std::cout << " std::unorderd_map size is taken into account \n";
    print_memory_usage();
  }

  std::cout << "================================================\n" ;
  std::cout << "after exiting the dictionary\n";
  print_memory_usage();
  std::cout << "================================================\n" ;
  
  { /* https://probablydance.com/2017/02/26/i-wrote-the-fastest-hashtable/ */
    ska::flat_hash_map<uint64_t, uint8_t> dictionary;

    timing = wtime();
    for (auto i = inp.begin(); i != inp.end(); ++i) {
      dictionary.emplace(std::make_pair(*i, 1));
    }
    timing = wtime() - timing;
    std::cout << "- insert: ska::unordered_map took " <<  timing << "sec \n";
    std::cout << "i.e. 2^" << std::log2(n/timing) << "elm/sec \n";

    uint64_t sum = 0;
    timing = wtime();
    for (auto i = qur.begin(); i != qur.end(); ++i) {
      auto found = dictionary.find(*i);
      if (found != dictionary.end()) 
	sum += found->first;
    }
    timing = wtime() - timing;
    std::cout << "- lookup: ska::flat_hash_map took " <<  timing << "sec \n";
    std::cout << "i.e. 2^" << std::log2(n/timing) << "elm/sec \n";
    std::cout << "- size = " << dictionary.max_size() << " bytes\n";
    std::cout << "sum = " << sum << "\n";
    std::cout << " ska::flat_hasp_map size is taken into account \n";
    print_memory_usage();

  }
  std::cout << "================================================\n" ;
  std::cout << "after exiting the dictionary\n";
  print_memory_usage();
  std::cout << "================================================\n" ;

}

