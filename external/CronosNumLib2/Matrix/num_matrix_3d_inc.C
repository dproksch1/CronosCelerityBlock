#include "matrix_3d_inc.C"

// specific 3D methods for NumMatrix
// -----------------------------------------------------------------------

template< >
__TYPE NumMatrix<__TYPE, 3>::
  interpol3(const int* pos, const int* diff, int r) const
{
  int i, j, k;
  __TYPE result = 0;

  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++)
      for (k = 0; k <= 3; k++)
	result += __TYPE(interpol3_w(i, double(diff[0])/r+1) 
	             * interpol3_w(j, double(diff[1])/r+1)
	             * interpol3_w(k, double(diff[2])/r+1) * 
          (*this)(pos[0]+i-1, pos[1]+j-1, pos[2]+k-1));

  return(result);
}

template<>
NumMatrix<__TYPE,3>& NumMatrix<__TYPE,3>::operator+=(const NumBoundary<__TYPE,3>& b)
{
  Boundary3d_forall(b, (*this)(x,y,z) += b(x,y,z);)
  return(*this);  
}

template<>
NumMatrix<__TYPE,3>& NumMatrix<__TYPE,3>::operator-=(const NumBoundary<__TYPE,3>& b)
{
  Boundary3d_forall(b, (*this)(x,y,z) -= b(x,y,z);)
  return(*this);  
}
