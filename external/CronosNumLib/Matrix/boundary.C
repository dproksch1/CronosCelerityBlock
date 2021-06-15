#include "boundary.H"

#define __TYPE int
#include "num_boundary_1d_inc.C"
#include "num_boundary_2d_inc.C"
#include "num_boundary_3d_inc.C"
#undef __TYPE
#define __TYPE float
#include "num_boundary_1d_inc.C"
#include "num_boundary_2d_inc.C"
#include "num_boundary_3d_inc.C"
#undef __TYPE
#define __TYPE double
#include "num_boundary_1d_inc.C"
#include "num_boundary_2d_inc.C"
#include "num_boundary_3d_inc.C"
#undef __TYPE
#define __TYPE double_complex
#include "num_boundary_1d_inc.C"
#include "num_boundary_2d_inc.C"
#include "num_boundary_3d_inc.C"
#undef __TYPE

template<class T, int rank>
NumBoundary<T, rank>::NumBoundary()
{
  index = 0;
}

template<class T, int rank>
NumBoundary<T, rank>::NumBoundary(const NumMatrix<T,rank>& matr, int width_)
{
  index = 0;

  newData(matr.getLow(),matr.getHigh(),width_);
}

template<class T, int rank> 
NumBoundary<T, rank>::~NumBoundary()
{
  deleteData();
}

template<class T, int rank> 
const int* NumBoundary<T, rank>::getLow() const
{
  return lo;
}

template<class T, int rank> 
const int* NumBoundary<T, rank>::getHigh() const
{
  return hi;
}

template<class T, int rank> 
const int* NumBoundary<T, rank>::getDims() const
{
  return dims;
}

template<class T, int rank>
int NumBoundary<T, rank>::getLow(int i) const
{
  return lo[i];
}

template<class T, int rank>
int NumBoundary<T, rank>::getHigh(int i) const
{
  return hi[i];
}

template<class T, int rank>
int NumBoundary<T, rank>::getDims(int i) const
{
  return dims[i];
}

template<class T, int rank>
int NumBoundary<T, rank>::getWidth() const
{
  return width;
}

template<class T, int rank>
void NumBoundary<T, rank>::resize(NumMatrix<T,rank>& matr, int width_)
{
  deleteData();
  newData(matr.getLow(),matr.getHigh(),width_);
}

// explicit instantiation

#define instNumBoundary(T, rank) \
template class NumBoundary<T, rank>

instNumBoundary(int, 1);
instNumBoundary(int, 2);
instNumBoundary(int, 3);
instNumBoundary(float, 1);
instNumBoundary(float, 2);
instNumBoundary(float, 3);
instNumBoundary(double, 1);
instNumBoundary(double, 2);
instNumBoundary(double, 3);
instNumBoundary(double_complex, 1);
instNumBoundary(double_complex, 2);
instNumBoundary(double_complex, 3);

#undef instNumBoundary



