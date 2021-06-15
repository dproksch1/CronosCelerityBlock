#include "matrix_1d_inc.H"

template< >
inline __TYPE NumMatrix<__TYPE,1>::mean(int i) const
{
  __TYPE res = __TYPE(0.25*((*this)(i+1)+(*this)(i-1))+0.5*(*this)(i));
  return res;
}

template< >
inline __TYPE NumMatrix<__TYPE,1>::mean(const int* pos) const
{
  return mean(pos[0]);
}

template< >
__TYPE NumMatrix<__TYPE, 1>::
interpol3(const int* pos, const int* diff, int r) const;

template<>
NumMatrix<__TYPE,1>& NumMatrix<__TYPE,1>::operator+=(const NumBoundary<__TYPE,1>& b);

template<>
NumMatrix<__TYPE,1>& NumMatrix<__TYPE,1>::operator-=(const NumBoundary<__TYPE,1>& b);
