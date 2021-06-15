//#ifndef FFTMATRIX_H
//#define FFTMATRIX_H
//
//#include "matrix.H"
//#include <rfftw.h>
//
//// Right now, this only works for
//// o  T == double
//// o  rank == 2
//
//
//template<class T, int rank>
//class FFTMatrix : public NumMatrix<T, rank> {
//public:
//  /** default constructur */
//  FFTMatrix();
//
//  /** constructor, which builds NumMatrix of size d[0] x ... x d[rank-1]
//      @memo constructor */
//  FFTMatrix(const int d[rank]);
//
//    /** constructor, which builds NumMatrix with lower 
//	indices l[0],...,l[rank-1]
//	and upper indices h[0],...,h[rank-1]
//	@memo constructor */
//  FFTMatrix(const int l[rank], const int h[rank]);
//
//  /** copy constructor */
//  FFTMatrix(const Matrix<T, rank>&);
//
//  /** index operator, writing */
//  T& real(const int* pos); // write
//  /** index operator, reading */
//  T  real(const int* pos) const; // read
//  /** index operator, writing */
//  T& real(int i);
//  /** index operator, reading */
//  T  real(int i) const;
//  /** index operator, writing */
//  T& real(int i, int j);
//  /** index operator, reading */
//  T  real(int i, int j) const;
//  /** index operator, writing */
//  T& real(int i, int j, int k);
//  /** index operator, reading */
//  T  real(int i, int j, int k) const;
//
//  /** index operator, writing */
//  double_complex& cmplx(const int* pos); // write
//  /** index operator, reading */
//  double_complex  cmplx(const int* pos) const; // read
//  /** index operator, writing */
//  double_complex& cmplx(int i);
//  /** index operator, reading */
//  double_complex  cmplx(int i) const;
//  /** index operator, writing */
//  double_complex& cmplx(int i, int j);
//  /** index operator, reading */
//  double_complex  cmplx(int i, int j) const;
//  /** index operator, writing */
//  double_complex& cmplx(int i, int j, int k);
//  /** index operator, reading */
//  double_complex  cmplx(int i, int j, int k) const;
//
//  void crFFT();
//  void rcFFT();
//
//protected:
//  virtual void newData(const int* l, const int* h);
//  virtual void deleteData();
//  
//private:
//   rfftwnd_plan rc_plan, cr_plan;
//};
//
//// -------------------------------------------------------------
//// inline functions
//
//template<class T, int rank>
//inline  T& FFTMatrix<T, rank>::real(int i)
//{
//  return NumMatrix<T, rank>::operator()(i);
//}
//
//template<class T, int rank>
//inline  T FFTMatrix<T, rank>::real(int i) const
//{
//  return NumMatrix<T, rank>::operator()(i);
//}
//
//template<class T, int rank>
//inline  T& FFTMatrix<T, rank>::real(int i, int j)
//{
//  return NumMatrix<T, rank>::operator()(i, j);
//}
//
//template<class T, int rank>
//inline  T FFTMatrix<T, rank>::real(int i, int j) const
//{
//  return NumMatrix<T, rank>::operator()(i, j);
//}
//
//template<class T, int rank>
//inline  T& FFTMatrix<T, rank>::real(int i, int j, int k)
//{
//  return NumMatrix<T, rank>::operator()(i, j, k);
//}
//
//template<class T, int rank>
//inline  T FFTMatrix<T, rank>::real(int i, int j, int k) const
//{
//  return NumMatrix<T, rank>::operator()(i, j, k);
//}
//
//template<class T, int rank>
//inline  double_complex& FFTMatrix<T, rank>::cmplx(int i)
//{
//  return (double_complex&)NumMatrix<T, rank>::operator()(2*i);
//}
//
//template<class T, int rank>
//inline  double_complex FFTMatrix<T, rank>::cmplx(int i) const
//{
//  return (double_complex)NumMatrix<T, rank>::operator()(2*i);
//}
//
//template<class T, int rank>
//inline  double_complex& FFTMatrix<T, rank>::cmplx(int i, int j)
//{
//  return (double_complex&)NumMatrix<T, rank>::operator()(2*i, j);
//}
//
//template<class T, int rank>
//inline  double_complex FFTMatrix<T, rank>::cmplx(int i, int j) const
//{
//  return (double_complex)NumMatrix<T, rank>::operator()(2*i, j);
//}
//
//template<class T, int rank>
//inline  double_complex& FFTMatrix<T, rank>::cmplx(int i, int j, int k)
//{
//  return (double_complex&)NumMatrix<T, rank>::operator()(2*i, j, k);
//}
//
//template<class T, int rank>
//inline  double_complex FFTMatrix<T, rank>::cmplx(int i, int j, int k) const
//{
//  return (double_complex)NumMatrix<T, rank>::operator()(2*i, j, k);
//}
//
//
//#endif // FFTMATRIX_H
