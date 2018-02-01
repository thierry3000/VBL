
#include <iostream>
#include <algorithm>
#include <vector>
#include "tbb/parallel_for_each.h"
#include "tbb/tick_count.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

int rnd() {
return (rand() % 100);
}

int main()
{
int addition = 5;
tbb::tick_count before = tbb::tick_count::now();
std::vector<int> numbers(100000000);
std::generate(numbers.begin(), numbers.end(), rnd);
tbb::tick_count after = tbb::tick_count::now();
std::printf("time spent (generate random): t%g secondsn", (after - before).seconds());

//std::for_each version
before = tbb::tick_count::now();
std::for_each(numbers.begin(), numbers.end(), [&] (int &v) { v += addition;});
after = tbb::tick_count::now();
std::printf("time spent (std::for_each): t%g secondsn", (after - before).seconds());

//tbb::parallel_for_each version
before = tbb::tick_count::now();
tbb::parallel_for_each(numbers.begin(), numbers.end(), [&] (int &v) { v += addition;});
after = tbb::tick_count::now();
std::printf("time spent (tbb::parallel_for_each): t%g secondsn", (after - before).seconds());

//tbb::parallel_for iterator version
before = tbb::tick_count::now();
tbb::parallel_for(
tbb::blocked_range<std::vector<int>::iterator>(numbers.begin(),
numbers.end()),
[&] (tbb::blocked_range<std::vector<int>::iterator> number) {
for (std::vector<int>::iterator it = number.begin(); it != number.end(); it++) {
(*it) += addition;
}
});
after = tbb::tick_count::now();
std::printf("time spent (tbb::parallel_for iterator): t%g secondsn", (after - before).seconds());

//tbb::parallel_for index version
before = tbb::tick_count::now();
tbb::parallel_for(size_t(0), size_t(numbers.size()),
[&] (size_t index) { numbers[index] += addition; });
after = tbb::tick_count::now();
std::printf("time spent (tbb::parallel_for index): t%g secondsn",(after - before).seconds());
printf("starting thierry test:\n");
return 0;
}
