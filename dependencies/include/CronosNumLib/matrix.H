#ifndef MATRIX_H
#define MATRIX_H

#include "complex.H"
#include "array.H"
#include <util.H>
#include <memory>
#include <array>

// ---------------------------------------------------------
// declarations

template<class T, int rank>
class Matrix;

template<class T, int rank>
std::ostream& operator<< (std::ostream&, const Matrix<T, rank>&);

/** An elementary matrix class */
template<class T, int rank>
class Matrix {
protected:
  T* matr;
  T* matr_fast;
  void* index;
  int size;
  std::array<int,rank> lo;
  std::array<int,rank> hi;
  std::array<int,rank> dims;
  std::string name;

public:
  /** default constructur */
  Matrix();

  /** constructor, which builds Matrix of size d[0] x ... x d[rank-1]

      Example:
      \begin{verbatim}
      int size[2]={512,512};
      Matrix<double,2> m(size);
      \end{verbatim}
      @memo constructor */
  Matrix(const int d[rank]);

  /** constructor, which builds Matrix with lower indices l[0],...,l[rank-1]
      and upper indices h[0],...,h[rank-1]

      Example:
      \begin{verbatim}
      int l[2]={-5,-10};
      int h[2]={15,36}
      Matrix<double,2> m(l,h);
      \end{verbatim}
      @memo constructor */
  Matrix(const int l[rank], const int h[rank]);

  /** copy constructor */
  Matrix(const Matrix<T, rank>&);

  /** move constructor */
  Matrix(Matrix<T, rank>&&);

  /** destructor */
  virtual ~Matrix();

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
  std::string getName() const;
  /** */
  void rename(const std::string &);
  /** */
  int getSize() const;
  /** return pointer to the first element holding the data,
      modifying allowed */
  T* data();
  /** return pointer to the first element holding the data,
      modifying forbidden */
  const T* data() const;
  /** return pointer to the first element holding the data,
      modifying allowed */
  operator T*();
  /** return pointer to the first element holding the data,
      modifying forbidden */
  operator const T*() const;

  /** resizes matrix to size d[0] x ... x d[rank-1] */
  void resize(const int* d);
  /** resizes to matrix with lower indices l[0],...,l[rank-1]
      and upper indices h[0],...,h[rank-1] */
  void resize(const int* l, const int* h);
  /** resizes to values of matrix m */
  void resize(const Matrix<T, rank>& m);

  /** index operator, writing */
  T& operator[](const int* pos); // write
  /** index operator, reading */
  T  operator[](const int* pos) const; // read
  /** index operator, writing */
  T& operator()(int i);
  /** index operator, reading */
  T  operator()(int i) const;
  /** index operator, writing */
  T& operator()(int i, int j);
  /** index operator, reading */
  T  operator()(int i, int j) const;
  /** index operator, writing */
  T& operator()(int i, int j, int k);
  /** index operator, reading */
  T  operator()(int i, int j, int k) const;
	/** index operator, writing */
	T& operator()(int i, int j, int k, int l);
	/** index operator, reading */
	T  operator()(int i, int j, int k, int l) const;
  /** index operator, writing */
  T& operator()(int i, int j, int k, int l, int m);
  /** index operator, reading */
  T  operator()(int i, int j, int k, int l, int m) const;

  /** assign */
  Matrix<T, rank>& operator=(const Matrix<T, rank>&);

  /** move assign */
  Matrix<T, rank>& operator=(Matrix<T, rank>&&);

  void swap(Matrix<T, rank>& other);

  /** comparison */
  int operator==(const Matrix<T, rank>&) const;
  /** comparison */
  int operator!=(const Matrix<T, rank>&) const;

  /** output to ostream */
  template<class T2, int rank2>
  friend std::ostream& operator<< (std::ostream&, const Matrix<T2, rank2>&);

protected:
  /** */
  virtual void deleteData();
  /** */
  virtual void newData(const int* l, const int* h);
};

// ---------------------------------------------------------

template<class T, int rank>
class NumBoundary;

/** An matrix class which can do
    elementary mathematical operations */
template<class T, int rank>
class NumMatrix : public Matrix<T, rank> {
private:
	int bcTypeLow[rank], bcTypeHigh[rank];
public:
  /** default constructur */
  NumMatrix();

  /** constructor, which builds NumMatrix of size d[0] x ... x d[rank-1]
      @memo constructor */
  NumMatrix(const int d[rank]);

    /** constructor, which builds NumMatrix with lower
	indices l[0],...,l[rank-1]
	and upper indices h[0],...,h[rank-1]
	@memo constructor */
  NumMatrix(const int l[rank], const int h[rank]);

  /** copy constructor */
  NumMatrix(const Matrix<T, rank>&);

  /** setting all elements of NumMatrix to zero */
  void clear();
  /** maximum norm of all elements of NumMatrix */
  double max_norm();
	
	// /** return max value of NumMatrx */
	// T get_max() const;
	// /** return min value of NumMatrx */
	// T get_min() const;

	/** return array holding lower boundary type*/
	const int* get_bcTypeLow() const;
	/** return array holding upper boundary type*/
	const int* get_bcTypeHigh() const;
	/** return lower boundary type by dimension */
	int get_bcTypeLow(int) const;
	/** return upper boundary type by dimension */
	int get_bcTypeHigh(int) const;
	/** set boundary type via array*/
	void set_bcType(const int*, const int*);
	/** set boundary type via array*/
	void set_bcType(int dir, int bcTL, int bcTH);
	/** set lower boundary type by dimension*/
	void set_bcTypeLow(int dir, int val);
	/** set upper boundary type by dimension*/
	void set_bcTypeHigh(int dir, int val);
	/** introduced by RK: - set constant value for whole field */
	void set_constVal(const T val);
	
  /** mean value of direct neighborhood */
  T mean(const int *pos) const;
  /** mean value of direct neighborhood */
  T mean(int i) const;
  /** mean value of direct neighborhood */
  T mean(int i, int j) const;
  /** mean value of direct neighborhood */
  T mean(int i, int j, int k) const;

  /** interpolates value at pos+diff/r (linear)
      @memo linear interpolation */
  T interpol(const int* pos, const int* diff, int r) const;
  /** interpolates value at pos+diff/r (cubic)
      @memo cubic interpolation */
  T interpol3(const int* pos, const int* diff, int r) const;

  /** */
  NumMatrix<T, rank>& operator += (const NumMatrix<T, rank>&);
  /** */
  NumMatrix<T, rank>& operator -= (const NumMatrix<T, rank>&);
  /** introduced by RK: */
  NumMatrix<T, rank>& operator *= (const NumMatrix<T, rank>&);
  /** introduced by RK: */
  NumMatrix<T, rank>& operator /= (const NumMatrix<T, rank>&);
  /** */
  NumMatrix<T, rank>& operator *= (T);
  /** */
  NumMatrix<T, rank>& operator /= (T);

  /** */
  NumMatrix<T, rank> operator+(const NumMatrix<T, rank>&) const;
  /** */
  NumMatrix<T, rank> operator-(const NumMatrix<T, rank>&) const;
  /** introduced by RK: */
//   NumMatrix<T, rank> operator*(const NumMatrix<T, rank>&) const;
  const NumMatrix<T, rank> operator*(const NumMatrix<T, rank>&) const;
  /** introduced by RK: */
  NumMatrix<T, rank> operator/(const NumMatrix<T, rank>&) const;
  /** */
  NumMatrix<T, rank> operator*(T) const;
  /** */
  NumMatrix<T, rank> operator/(T) const;

  /** */
  NumMatrix<T, rank>& operator+=(const NumBoundary<T, rank>&);
  /** */
  NumMatrix<T, rank>& operator-=(const NumBoundary<T, rank>&);
};

// -------------------------------------------------------------
// inline functions

// write function
template<class T, int rank>
inline T& Matrix<T, rank>::operator[](const int* pos)
{
  assert(pos[rank-1] >= lo[rank-1]);
  assert(pos[rank-1] <= hi[rank-1]);
  // std::cout << " Ich schreibe " << pos[0] << std::endl;
  int p = pos[rank-1]-lo[rank-1];


  for(int i = rank-2; i >= 0; i--) {
    assert(pos[i] >= lo[i]) ;
    assert(pos[i] <= hi[i]);
    p = p*dims[i] + (pos[i]-lo[i]);
  }

  assert((p >= 0) && (p < size));
  return matr[p];
}

// read function
template<class T, int rank>
inline T Matrix<T, rank>::operator[](const int* pos) const
{
  assert(pos[rank-1] >= lo[rank-1]);
  assert(pos[rank-1] <= hi[rank-1]);
  // std::cout << " Here I am " << pos[0] << std::endl;
  int p = pos[rank-1];

  // if(rank==1 && pos[0]==32) {
  // }

  for(int i = rank-2; i >= 0; i--) {
    assert(pos[i] >= lo[i]);
    assert(pos[i] <= hi[i]);
    p = p*dims[i] + pos[i];
  }

  return matr_fast[p];
}

template<class T, int rank>
inline T* Matrix<T, rank>::data()
{
  return matr;
}

template<class T, int rank>
inline const T* Matrix<T, rank>::data() const
{
  return matr;
}

template<class T, int rank>
inline Matrix<T, rank>::operator T*()
{
  return data();
}

template<class T, int rank>
inline Matrix<T, rank>::operator const T*() const
{
  return data();
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i)
{
  assert(rank == 1);
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  // if(i>26 && i<36)
  // 	std::cout << " Andere Runde Klammer " << i << std::endl;
  return matr_fast[i];
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i) const
{
	// if(i>26 && i<36)
	// 	std::cout << " Runde Klammer " << i << std::endl;
  assert(rank == 1);
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  return matr_fast[i];
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i, int j)
{
  assert(rank == 2);
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  return matr_fast[i + dims[0] * j];
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j) const
{
  assert(rank == 2);
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  std::cerr << "()" << i+dims[0]*j << std::endl;
  return matr_fast[i + dims[0] * j];
}

template<class T, int rank>
inline  T& Matrix<T, rank>::operator ()(int i, int j, int k)
{
  assert(rank == 3);
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  assert(k >= lo[2]);
  assert(k <= hi[2]);
  return matr_fast[i + dims[0]*(j + dims[1]*k)];
}

template<class T, int rank>
inline  T Matrix<T, rank>::operator ()(int i, int j, int k) const
{
  assert(rank == 3);
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  assert(k >= lo[2]);
  assert(k <= hi[2]);
  return matr_fast[i + dims[0]*(j + dims[1]*k)];
}

template<class T, int rank>
inline void Matrix<T, rank>::swap(Matrix<T, rank>& other)
{
  using std::swap;
  swap(this->matr, other.matr);
  swap(this->matr_fast, other.matr_fast);
  swap(this->index, other.index);
  swap(this->size, other.size);
  swap(this->lo, other.lo);
  swap(this->hi, other.hi);
  swap(this->dims, other.dims);
  swap(this->name, other.name);
}

template<class T, int rank>
inline Matrix<T, rank>& Matrix<T, rank>::operator=(Matrix<T, rank>&& other)
{
  Matrix<T, rank> tmp(other);
  this->swap(tmp);
  return *this;
}

template<class T, int rank>
inline Matrix<T, rank>::Matrix(Matrix<T, rank>&& other): Matrix()
{
  other.swap(*this);
}


// Matrix Conversion

template<class To, class From, int rank>
Matrix<To, rank> convertMatrix(const Matrix<From, rank> & matr_in) {
  Matrix<To, rank> matr_out(matr_in.getLow(), matr_in.getHigh());
  for (int i = 0; i < matr_in.getSize(); i++)
    matr_out.data()[i] = static_cast<To>(matr_in.data()[i]);

  return matr_out;
}




#define __TYPE char
#include "matrix_1d_inc.H"
#include "matrix_2d_inc.H"
#include "matrix_3d_inc.H"
#include "matrix_4d_inc.H"
#include "matrix_5d_inc.H"
#undef __TYPE
#define __TYPE int
#include "num_matrix_1d_inc.H"
#include "num_matrix_2d_inc.H"
#include "num_matrix_3d_inc.H"
#include "num_matrix_4d_inc.H"
#include "num_matrix_5d_inc.H"
#undef __TYPE
#define __TYPE float
#include "num_matrix_1d_inc.H"
#include "num_matrix_2d_inc.H"
#include "num_matrix_3d_inc.H"
#include "num_matrix_4d_inc.H"
#include "num_matrix_5d_inc.H"
#undef __TYPE
#define __TYPE double
#include "num_matrix_1d_inc.H"
#include "num_matrix_2d_inc.H"
#include "num_matrix_3d_inc.H"
#include "num_matrix_4d_inc.H"
#include "num_matrix_5d_inc.H"
#undef __TYPE
#define __TYPE double_complex
#include "num_matrix_1d_inc.H"
#include "num_matrix_2d_inc.H"
#include "num_matrix_3d_inc.H"
#include "num_matrix_4d_inc.H"
#include "num_matrix_5d_inc.H"
#undef __TYPE


#endif
