#include "matrix_1d_inc.C"

// specific 1D methods for NumMatrix
// -----------------------------------------------------------------------

template< >
__TYPE NumMatrix<__TYPE, 1>::
  interpol3(const int* pos, const int* diff, int r) const
{
  int i;
  __TYPE result = 0;

  for (i = 0; i <= 3; i++)
      result += __TYPE(interpol3_w(i, double(diff[0])/r+1)
	* (*this)(pos[0]+i-1));

  return(result);
} 

template<>
NumMatrix<__TYPE,1>& NumMatrix<__TYPE,1>::operator+=(const NumBoundary<__TYPE,1>& b)
{
  Boundary1d_forall(b, (*this)(x) += b(x);)
  return(*this);  
}

template<>
NumMatrix<__TYPE,1>& NumMatrix<__TYPE,1>::operator-=(const NumBoundary<__TYPE,1>& b)
{
  Boundary1d_forall(b, (*this)(x) -= b(x);)
  return(*this);  
}
