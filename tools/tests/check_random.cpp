/* I encoutered problems with random numbers on different mashines,
 * so I try to test the environment with this.
 */
# include <boost/random.hpp>
# include <numeric> // std::accumulate
boost::random::mt19937 rng;         // produces randomness out of thin air
                                    // see pseudo-random number generators
boost::random::uniform_01<> six_f;

int main(int argc, char *argv[])
{
  //works for int, see boost examples
  
  
  boost::random::uniform_int_distribution<> six(1,6);
                                    // distribution that maps to 1..6
                                    // see random number distributions
  for( int i=0;i<100;i++)
  {
    int x = six(rng);                   // simulate rolling a die
    printf("i: %i, %i\n", i,x);
  }
  
//   boost::random::uniform_real_distribution<> six_f(1,6);
  
  const int n=10000;
//   for( int i=0;i<100;i++)
//   {
//     float x = six_f(rng);                   // simulate rolling a die
//     printf("i: %i, %0.2f\n", i,x);
//   }
  std::vector<float> my_random_numbers;
  for( int i=0;i<n;i++)
  {
    float x = six_f(rng);                   // simulate rolling a die
    my_random_numbers.push_back(x);
  }
  double sum = std::accumulate(my_random_numbers.begin(), my_random_numbers.end(), 0.0);
  double mean = sum / my_random_numbers.size();

  double sq_sum = std::inner_product(my_random_numbers.begin(), my_random_numbers.end(), my_random_numbers.begin(), 0.0);
  double var = sq_sum / my_random_numbers.size() - mean * mean;
  printf("average: %0.2f (should 0.5) \n", mean);
  printf("var: %0.2f (should %0.2f) \n", var, (1./sqrt(n)));
}
