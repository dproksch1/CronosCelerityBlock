#include "matrix_2d_inc.C"

// specific 2D methods for NumMatrix
// -----------------------------------------------------------------------

template< >
__TYPE NumMatrix<__TYPE, 2>::
  interpol3(const int* pos, const int* diff, int r) const
{
  int i, j;
  __TYPE result = 0;

  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++)
      result += __TYPE(interpol3_w(i, double(diff[0])/r+1) 
  	           * interpol3_w(j, double(diff[1])/r+1) 
	* (*this)(pos[0]+i-1, pos[1]+j-1));

  return(result);
}

template<>
NumMatrix<__TYPE,2>& NumMatrix<__TYPE,2>::operator+=(const NumBoundary<__TYPE,2>& b)
{
  Boundary2d_forall(b, (*this)(x,y) += b(x,y);)
  return(*this);  
}

template<>
NumMatrix<__TYPE,2>& NumMatrix<__TYPE,2>::operator-=(const NumBoundary<__TYPE,2>& b)
{
  Boundary2d_forall(b, (*this)(x,y) -= b(x,y);)
  return(*this);  
}

// // #if (__TYPE != double_complex)
// template<>
// __TYPE NumMatrix<__TYPE,2>::get_min() const
// {
// 	__TYPE minimum = this->matr[0];
// 	for (int i = 0; i < this->size; i++) {
// 		minimum = std::min(minimum, this->matr[i]);
// 	}
// 	return minimum;
// }

// template< >
// inline __TYPE NumMatrix<__TYPE,2>::get_max() const
// {
// 	__TYPE maximum = this->matr[0];
// 	for (int i = 0; i < this->size; i++) {
// 		maximum = std::max(maximum, this->matr[i]);
// 	}
// 	return maximum;
// }
// // #endif
