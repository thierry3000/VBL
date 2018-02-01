#include "tbb/tbb.h"
#include "tbb/parallel_do.h"

#include <list>
#include <chrono>

#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
boost::random::random_device gen;
boost::random::uniform_01<> uni_f;

double myRandom()
{
  return uni_f(gen);
}

void Foo(const float aFloat)
{
  //printf(" doing a complicated new calculation on: %f\n", aFloat);
  const int zahl = 8e4 ;
  std::vector<double> vec_of_rand;
  int i=0;
  for (int p=0;p<zahl;p++)
  {
    double ran = myRandom();
    vec_of_rand.push_back(ran);
  }
  std::sort(vec_of_rand.begin(), vec_of_rand.end());
}
/*********************
 * SERIAL WAY
 **************************/
void SerialApplyFoo( float a[], size_t n ) {
    for( size_t i=0; i!=n; ++i )
        Foo(a[i]);
}

/*********************
 * CLASS 
 * using tbb::blocked_range
 * this is more or less openmp
 **************************/
class ApplyFoo {
    float *const my_a;
public:
    void operator()( const tbb::blocked_range<size_t>& r ) const {
        float *a = my_a;
        for( size_t i=r.begin(); i!=r.end(); ++i ) 
           Foo(a[i]);
    }
    ApplyFoo( float a[] ) :
        my_a(a)
    {}
};

void ParallelApplyFoo( float a[], size_t n ) {
    parallel_for(tbb::blocked_range<size_t>(0,n), ApplyFoo(a));
}

/*********************
 * CLASS 
 * using tbb::parallel_do
 * this is works for unknow number of entries
 * ---> run until done!
 **************************/
class ApplyFooList{
public:
  ApplyFooList(){};
  void operator()( const float& item) const
  {
    Foo(item);
  };
};

void SerialApplyFooToList( const std::list<float>& list)
{
  for( std::list<float>::const_iterator i=list.begin(); i!=list.end(); ++i)
  {
    Foo(*i);
  }
}
void ParallelApplyFooToList( const std::list<float> &list)
{
  tbb::parallel_do(list.begin(), list.end(), ApplyFooList());
}
int main()
{
  const size_t n=10;
  const size_t n2=100;
  float bar[10] = {42,43,44,45,43,43,43,45,656,7};
  std::vector<float> myRands;
  std::list<float> myRands2;
  for(size_t i = 0; i<n2; i++)
  {
    myRands.push_back(myRandom());
    myRands2.push_back(myRandom());
  }
  auto start = std::chrono::high_resolution_clock::now();
  SerialApplyFoo(&myRands[0], myRands.size());
  auto stop = std::chrono::high_resolution_clock::now();
  auto SerialApplyFoo = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  start = std::chrono::high_resolution_clock::now();
  ParallelApplyFoo(&myRands[0], myRands.size());
  stop = std::chrono::high_resolution_clock::now();
  auto ParallelApplyFoo = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  
  start = std::chrono::high_resolution_clock::now();
  SerialApplyFooToList(myRands2);
  stop = std::chrono::high_resolution_clock::now();
  auto SerialApplyFooToListDuration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  
  start = std::chrono::high_resolution_clock::now();
  ParallelApplyFooToList(myRands2);
  stop = std::chrono::high_resolution_clock::now();
  auto ParallelApplyFooToListDuration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  
  std::cout << "Time taken serial  : "
         << SerialApplyFoo.count() << " microseconds" << std::endl;
  std::cout << "Time taken parallel: "
         << ParallelApplyFoo.count() << " microseconds" << std::endl;
         
  std::cout << "Time taken serial  : "
         << SerialApplyFooToListDuration.count() << " microseconds" << std::endl;
  std::cout << "Time taken parallel: "
         << ParallelApplyFooToListDuration.count() << " microseconds" << std::endl;
//   for(size_t i = 0; i<n; i++)
//   {
//     printf("i: %i, bar[i]: %f\n" , i, bar[i]);
//   }
//   SerialApplyFoo(bar, n);
}

