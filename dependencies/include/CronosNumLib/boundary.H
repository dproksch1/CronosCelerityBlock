#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "matrix.H"

/** Boundary class */
template<class T, int rank>
class NumBoundary {
  void* index;
  int width;
  int lo[rank];
  int hi[rank];
  int dims[rank];

public:
  /** constructor */
  NumBoundary();
  /** copy constructor */
  NumBoundary(const NumBoundary<T,rank>&);
  /** constructs Boundary with size corresponding to m and width w */
  NumBoundary(const NumMatrix<T,rank>& m, int w);
  /** destructor */
  ~NumBoundary();

  /** */
  const int* getLow() const;
  /** */
  const int* getHigh() const;
  /** */
  const int* getDims() const;
  /** */
  int getLow(int) const;
  /** */
  int getHigh(int) const;
  /** */
  int getDims(int) const;
  /** */
  int getWidth() const;

  /** */
  void resize(NumMatrix<T,rank>&, int);
  /** */
  void clear();

  /** index operator, modifying allowed */
  T& operator[](const int* pos); // write
  /** index operator, modifying forbidden */
  T  operator[](const int* pos) const; // read
  /** index operator, modifying allowed */
  T& operator()(int i);
  /** index operator, modifying forbidden */
  T  operator()(int i) const;
  /** index operator, modifying allowed */
  T& operator()(int i, int j);
  /** index operator, modifying forbidden */
  T  operator()(int i, int j) const;
  /** index operator, modifying allowed */
  T& operator()(int i, int j, int k);
  /** index operator, modifying forbidden */
  T  operator()(int i, int j, int k) const;

  /** assign */
  NumBoundary& operator=(const NumBoundary<T,rank>&);
  /** */
  NumBoundary& operator+=(const NumMatrix<T,rank>&);
  /** */
  NumBoundary& operator-=(const NumMatrix<T,rank>&);
  /** */
  NumBoundary& operator/=(T);
  /** */
  NumBoundary& operator*=(T);
   
protected:
  /** */
  void newData(const int *l, const int *h, int width);
  /** */
  void deleteData();
};

// ------------------------------------------------------------
// inline functions

#define __TYPE int
#include "num_boundary_1d_inc.H"
#include "num_boundary_2d_inc.H"
#include "num_boundary_3d_inc.H"
#undef __TYPE
#define __TYPE float
#include "num_boundary_1d_inc.H"
#include "num_boundary_2d_inc.H"
#include "num_boundary_3d_inc.H"
#undef __TYPE
#define __TYPE double
#include "num_boundary_1d_inc.H"
#include "num_boundary_2d_inc.H"
#include "num_boundary_3d_inc.H"
#undef __TYPE
#define __TYPE double_complex
#include "num_boundary_1d_inc.H"
#include "num_boundary_2d_inc.H"
#include "num_boundary_3d_inc.H"
#undef __TYPE


#endif
