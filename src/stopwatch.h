#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_

#include <ctime>
#include <cstdio>
#include <string>
using std::string;

/** 
 * \file stopwatch.h
 * \brief Stopwatch class
 * \author D. M. Luchtenburg
 *
 * This class describes the implementation of the stopwatch
 */

class Stopwatch {
 public:
  Stopwatch();
  ~Stopwatch();
  void toc();
 private:
  clock_t start;
};

#endif // _STOPWATCH_H_
