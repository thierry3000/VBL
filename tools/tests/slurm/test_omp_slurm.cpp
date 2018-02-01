#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include <chrono>
#include <thread>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
boost::random::random_device gen;
boost::random::uniform_01<> uni_f;

double myRandom()
{
  return uni_f(gen);
}

int main (int argc, char *argv[]) 
{
int nthreads, tid;

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel private(nthreads, tid)
  {
  const size_t n2=1000000;
  std::vector<float> myRands2;
  for(size_t i = 0; i<n2; i++)
  {
    myRands2.push_back(myRandom());
  }
  std::sort(myRands2.begin(), myRands2.end());
  /* Obtain thread number */
  tid = omp_get_thread_num();
  printf("Hello World from thread = %d\n", tid);
  
  /* Only master thread does this */
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }

  }  /* All threads join master thread and disband */

}

