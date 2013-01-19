#include "stopwatch.h"
using std::string;

Stopwatch::Stopwatch() : start(std::clock()) {}

Stopwatch::~Stopwatch() {}

void Stopwatch::toc() {
  clock_t total = clock() - start;
  printf("Elapsed time is %.6e seconds.\n", double(total) / CLOCKS_PER_SEC);
}
