#include "matrix.H"
#include "boundary.H"

//#include <stl_algobase.h>
#include <algorithm>
#include <math.h>
#include <stdlib.h>

using namespace std;

// --------------------------------------------------------------
// implementation

template<class T, int rank>
Matrix<T, rank>::Matrix()
{
  matr = 0;
  size = 0;
  index = 0;
  for (int dir = 0; dir < rank; dir++) {
	  lo[dir] = 0;
	  hi[dir] = 0;
  }
  name = "void";
}

template<class T, int rank>
Matrix<T, rank>::Matrix(const int d[rank])
{
  matr = 0;
  size = 0;

  Index l(rank);
  Index h(rank);
  l.clear();
  for (int i = 0; i < rank; i++)
    h[i] = d[i] - 1;

  newData(l, h);
  name = "void";
}

template<class T, int rank>
Matrix<T, rank>::Matrix(const int l[rank], const int h[rank])
{
  matr = 0;
  size = 0;
  newData(l, h);
  name = "void";
}

template<class T, int rank>
Matrix<T, rank>::Matrix(const Matrix<T, rank>& _matr)
{
  matr = 0;
  size = 0;

  newData(_matr.lo.data(), _matr.hi.data());

  for (int i = 0; i < size; i++)
    matr[i] = _matr.matr[i];
  name = "void";
}

template<class T, int rank>
Matrix<T, rank>::~Matrix()
{
  deleteData();
}

template<class T, int rank>
const int* Matrix<T, rank>::getLow() const
{
  return lo.data();
}

template<class T, int rank>
const int* Matrix<T, rank>::getHigh() const
{
  return hi.data();
}

template<class T, int rank>
const int* Matrix<T, rank>::getDims() const
{
  return dims.data();
}

template<class T, int rank>
int Matrix<T, rank>::getLow(int i) const
{
  return lo[i];
}

template<class T, int rank>
int Matrix<T, rank>::getHigh(int i) const
{
  return hi[i];
}

template<class T, int rank>
int Matrix<T, rank>::getDims(int i) const
{
  return dims[i];
}

template<class T, int rank>
int Matrix<T, rank>::getSize() const
{
  return size;
}

template<class T, int rank>
string Matrix<T, rank>::getName() const
{
  return name;
}

template<class T, int rank>
void Matrix<T, rank>::rename(const string &name)
{
  this->name = name;
}

template<class T, int rank>
Matrix<T, rank>& Matrix<T, rank>::operator=(const Matrix<T, rank>& _matr)
{
  for (int d = 0; d < rank; d++) {
    if (_matr.lo[d] != lo[d] ||
	_matr.hi[d] != hi[d]) {
      deleteData();
      newData(_matr.lo.data(), _matr.hi.data());
      break;
    }
  }

  for (int i = 0; i < size; i++)
    matr[i] = _matr.matr[i];
  return *this;
}

template<class T, int rank>
int Matrix<T, rank>::operator==(const Matrix<T, rank>& _matr) const
{
  int i;
  for (i = 0; i < rank; i++) {
    if (lo[i] != _matr.lo[i])
      return 0;
    if (hi[i] != _matr.hi[i])
      return 0;
  }

  for (i = 0; i < size; i++)
    if (matr[i] != _matr.matr[i])
      return 0;
  return 1;
}

template<class T, int rank>
int Matrix<T, rank>::operator!=(const Matrix<T, rank>& _matr) const
{
  return !operator==(_matr);
}

template<class T, int rank>
void Matrix<T, rank>::resize(const int* d)
{
  Index l(rank);
  Index h(rank);
  l.clear();
  for (int i = 0; i < rank; i++)
    h[i] = d[i] - 1;

  resize(l, h);
}

template<class T, int rank>
void Matrix<T, rank>::resize(const int* l, const int* h)
{
  bool sameExtent = true;
  for (int dim = 0; dim < rank; ++dim) {
    if (l[dim] != lo[dim] || h[dim] != hi[dim]) {
      sameExtent = false;
      break;
    }
  }
  if (!sameExtent || size == 0) {
    deleteData();
    newData(l, h);
  }
}

template<class T, int rank>
void Matrix<T, rank>::resize(const Matrix<T, rank>& _matr)
{
  resize(_matr.lo.data(), _matr.hi.data());
}

template<class T, int rank>
void Matrix<T, rank>::deleteData()
{
  if (matr)
    delete[] matr;
  matr = 0;
  size = 0;
}

template<class T, int rank>
void Matrix<T, rank>::newData(const int* l, const int* h)
{
  size = 1;
  int d;
  for (d = 0; d < rank; d++) {
    assert(l[d] <= h[d]);
    lo[d] = l[d];
    hi[d] = h[d];
    dims[d] = h[d] - l[d] + 1;
    size *= dims[d];
  }
  matr = new T[size];
  int p = -lo[rank-1];

  for (d = rank-2; d >= 0 ; d--) {
    p = p*dims[d] -lo[d];
  }
  matr_fast = matr + p;
}

template<class T, int rank>
ostream& operator<<(ostream& os, const Matrix<T, rank>& matr)
{
  if (!matr) {
    os << "(undef)" << endl;
    return os;
  }
  os << "ostream& operator<<(ostream& os, const Matrix<T, rank>& matr)"
     << endl << " not implemented" << endl;
  return os;
}

// -----------------------------------------------------------------

template<class T, int rank>
NumMatrix<T, rank>::NumMatrix() : Matrix<T, rank>() {}

template<class T, int rank>
NumMatrix<T, rank>::NumMatrix(const int d[rank]) : Matrix<T, rank>(d) {}

template<class T, int rank>
NumMatrix<T, rank>::NumMatrix(const int l[rank], const int h[rank])
  : Matrix<T, rank>(l, h) {}

template<class T, int rank>
NumMatrix<T, rank>::NumMatrix(const Matrix<T, rank>& _matr)
  : Matrix<T, rank>(_matr) {}

template<class T, int rank>
void NumMatrix<T, rank>::clear()
{
  for (int i = 0; i < this->size; i++)
    this->matr[i] = 0;
}

template<class T, int rank>
double NumMatrix<T, rank>::max_norm()
{
  double maximum = abs(this->matr[0]);
  for (int i = 0; i < this->size; i++)
    maximum=std::max(maximum,double(abs(this->matr[i])));
  return maximum;
}


template<class T, int rank>
const int* NumMatrix<T, rank>::get_bcTypeLow() const
{
  return bcTypeLow;
}


template<class T, int rank>
const int* NumMatrix<T, rank>::get_bcTypeHigh() const
{
  return bcTypeHigh;
}


template<class T, int rank>
int NumMatrix<T, rank>::get_bcTypeLow(int dir) const
{
  return bcTypeLow[dir];
}


template<class T, int rank>
int NumMatrix<T, rank>::get_bcTypeHigh(int dir) const
{
  return bcTypeHigh[dir];
}


template<class T, int rank>
void NumMatrix<T, rank>::set_bcType(const int* bcTypeLow, const int* bcTypeHigh)
{
	for (int dir = 0; dir < rank; dir++) {
		this->bcTypeLow[dir] = bcTypeLow[dir];
		this->bcTypeHigh[dir] = bcTypeHigh[dir];
	}
}


template<class T, int rank>
void NumMatrix<T, rank>::set_bcType(int dir, int bcTypeLow, int bcTypeHigh)
{
	this->bcTypeLow[dir] = bcTypeLow;
	this->bcTypeHigh[dir] = bcTypeHigh;
}


template<class T, int rank>
void NumMatrix<T, rank>::set_bcTypeLow(int dir, int bcTypeLow)
{
	this->bcTypeLow[dir] = bcTypeLow;
}


template<class T, int rank>
void NumMatrix<T, rank>::set_bcTypeHigh(int dir, int bcTypeHigh)
{
	this->bcTypeHigh[dir] = bcTypeHigh;
}

template<class T, int rank>
void NumMatrix<T, rank>::set_constVal(const T val)
{
  for (int i = 0; i < this->size; i++)
    this->matr[i] = val;
}



template<class T, int rank>
T NumMatrix<T, rank>::interpol(const int* pos, const int* diff, int r) const
{
  int dual, i;
  Index p(rank);

  T res = 0;

  for (dual = 0; dual < (1 << rank); dual++) {
    int fak = 1;
    for (i = 0; i < rank; i++) {
      if (dual & (1 << i)) {
	fak *= r - diff[i];
	p[i] = pos[i];
      } else {
	fak *= diff[i];
	p[i] = pos[i] + 1;
      }
    }
    res += T(double(fak) * (*this)[p]);
  }

  for (i = 0; i < rank; i++)
    res /= r;
  return res;
}

// template<class T, int rank>
// T NumMatrix<T, rank>::interpol3(const int* pos, const int* diff, int r) const
// {
//   cerr << "fatal: interpol3 not generally implemented!" << endl;
//   return (T)0;
// }

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator*=(T t)
{
  for (int i = 0; i < this->size; i++)
    this->matr[i] *= t;
  return *this;
}

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator/=(T t)
{
  T div = T(1./t);
  for (int i = 0; i < this->size; i++)
    this->matr[i] *= div;
  return *this;
}

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator+=(const NumMatrix<T, rank>&
						   _matr)
{
  assert(this->size == _matr.size);
  for (int i = 0; i < this->size; i++)
    this->matr[i] += _matr.matr[i];
  return *this;
}

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator-=(const NumMatrix<T, rank>&
						   _matr)
{
  assert(this->size == _matr.size);
  for (int i = 0; i < this->size; i++)
    this->matr[i] -= _matr.matr[i];
  return *this;
}

/*Begin new: added by RK*/
template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator*=(const NumMatrix<T, rank>&
						   _matr)
{
  assert(this->size == _matr.size);
  for (int i = 0; i < this->size; i++)
    this->matr[i] *= _matr.matr[i];
  return *this;
}

template<class T, int rank>
NumMatrix<T, rank>& NumMatrix<T, rank>::operator/=(const NumMatrix<T, rank>&
						   _matr)
{
  assert(this->size == _matr.size);
  for (int i = 0; i < this->size; i++)
    this->matr[i] /= _matr.matr[i];
  return *this;
}
/*End new*/

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::
  operator+(const NumMatrix<T, rank>& matr1) const
{
  assert(this->size == matr1.size);
  NumMatrix<T, rank> tmp = (*this);
  tmp += matr1;
  return tmp;
}

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::
  operator-(const NumMatrix<T, rank>& matr1) const
{
  assert(this->size == matr1.size);
  NumMatrix<T, rank> tmp = (*this);
  tmp -= matr1;
  return tmp;
}

// /*Begin new: added by RK*/
// template<class T, int rank>
// NumMatrix<T, rank> NumMatrix<T, rank>::
//   operator*(const NumMatrix<T, rank>& matr1) const
// {
//   assert(size == matr1.size);
//   NumMatrix<T, rank> tmp = (*this);
//   tmp *= matr1;
//   return tmp;
// }
/*Begin new: added by RK*/
template<class T, int rank>
const NumMatrix<T, rank> NumMatrix<T, rank>::
  operator*(const NumMatrix<T, rank>& matr1) const
{
  assert(this->size == matr1.size);
  NumMatrix<T, rank> tmp = (*this);
  tmp *= matr1;
  return tmp;
}

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::
  operator/(const NumMatrix<T, rank>& matr1) const
{
  assert(this->size == matr1.size);
  NumMatrix<T, rank> tmp = (*this);
  tmp /= matr1;
  return tmp;
}
/*End new*/

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::operator*(T t) const
{
  NumMatrix<T, rank> tmp = (*this);
  tmp *= t;
  return tmp;
}

template<class T, int rank>
NumMatrix<T, rank> NumMatrix<T, rank>::operator/(T t) const
{
  NumMatrix<T, rank> tmp = (*this);
  tmp /= t;
  return tmp;
}


// ------------------------------------------------------------
// specializations

inline double interpol3_w(int i, double x)
{
  switch (i) {
  case 0: return((x*x*x-6*x*x+11*x-6)/(-6));
  case 1: return((x*x*x-5*x*x+6*x)/2);
  case 2: return((x*x*x-4*x*x+3*x)/(-2));
  case 3: return((x*x*x-3*x*x+2*x)/6);
  }
  return -1;
}

#define __TYPE char
#include "matrix_1d_inc.C"
#include "matrix_2d_inc.C"
#include "matrix_3d_inc.C"
#undef __TYPE
#define __TYPE int
#include "num_matrix_1d_inc.C"
#include "num_matrix_2d_inc.C"
#include "num_matrix_3d_inc.C"
#undef __TYPE
#define __TYPE float
#include "num_matrix_1d_inc.C"
#include "num_matrix_2d_inc.C"
#include "num_matrix_3d_inc.C"
#undef __TYPE
#define __TYPE double
#include "num_matrix_1d_inc.C"
#include "num_matrix_2d_inc.C"
#include "num_matrix_3d_inc.C"
#undef __TYPE
#define __TYPE double_complex
#include "num_matrix_1d_inc.C"
#include "num_matrix_2d_inc.C"
#include "num_matrix_3d_inc.C"
#undef __TYPE


// explicit instantiation

#define instMatrix(T, rank) \
template class Matrix<T, rank>; \
//template ostream& ::operator<<(ostream& os, const Matrix<T, rank>& matr)

#define instNumMatrix(T, rank) \
template class Matrix<T, rank>; \
template class NumMatrix<T, rank>; \
//template ostream& ::operator<<(ostream& os, const Matrix<T, rank>& matr)

instMatrix(char, 1);
instMatrix(char, 2);
instMatrix(char, 3);

instNumMatrix(bool, 1);
instNumMatrix(bool, 2);
instNumMatrix(bool, 3);
instNumMatrix(bool, 4);
instNumMatrix(int, 1);
instNumMatrix(int, 2);
instNumMatrix(int, 3);
instNumMatrix(int, 4);
instNumMatrix(long, 1);
instNumMatrix(long, 2);
instNumMatrix(long, 3);
instNumMatrix(long, 4);
instNumMatrix(float, 1);
instNumMatrix(float, 2);
instNumMatrix(float, 3);
instNumMatrix(float, 4);
instNumMatrix(double, 1);
instNumMatrix(double, 2);
instNumMatrix(double, 3);
instNumMatrix(double, 4);
instNumMatrix(double, 5);
instNumMatrix(double_complex, 1);
instNumMatrix(double_complex, 2);
instNumMatrix(double_complex, 3);

#undef instMatrix
#undef instNumMatrix

