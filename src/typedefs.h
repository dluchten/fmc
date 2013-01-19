#ifndef TYPEDEFS_H_
#define TYPEDEFS_H_

#include <vector>

// Data containers
typedef std::vector<double> VecDoub1D;
typedef std::vector<VecDoub1D> VecDoub2D;
typedef std::vector<VecDoub2D> VecDoub3D;
typedef std::vector<VecDoub3D> VecDoub4D;

typedef double (*DoubFuncInts)(const int *x);

#endif  // TYPEDEFS_H_
