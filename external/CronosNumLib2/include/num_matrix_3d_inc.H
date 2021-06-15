#include "matrix_3d_inc.H"

template< >
inline __TYPE NumMatrix<__TYPE,3>::mean(int i, int j, int k) const
{
  __TYPE res = __TYPE(1./12*((*this)(i+1,j,k)+(*this)(i-1,j,k)
               +(*this)(i,j+1,k)+(*this)(i,j-1,k)
	       +(*this)(i,j,k+1)+(*this)(i,j,k-1))
		  +0.5*(*this)(i,j,k));

//      __TYPE res = __TYPE(0.125*(*this)(i,j,k)
//            +0.0625*((*this)(i+1,j,k)+(*this)(i-1,j,k)
//                    +(*this)(i,j+1,k)+(*this)(i,j-1,k)
//     	           +(*this)(i,j,k+1)+(*this)(i,j,k-1))
//           +0.03125*((*this)(i+1,j+1,k)+(*this)(i-1,j+1,k)
// 		   +(*this)(i+1,j-1,k)+(*this)(i-1,j-1,k)
// 		   +(*this)(i,j+1,k+1)+(*this)(i,j-1,k+1)
// 		   +(*this)(i,j+1,k-1)+(*this)(i,j-1,k-1)
// 		   +(*this)(i+1,j,k+1)+(*this)(i+1,j,k-1)
// 		   +(*this)(i-1,j,k+1)+(*this)(i-1,j,k-1))
//          +0.015625*((*this)(i+1,j+1,k+1)+(*this)(i-1,j+1,k+1)
// 		   +(*this)(i+1,j-1,k+1)+(*this)(i-1,j-1,k+1)
// 		   +(*this)(i+1,j-1,k-1)+(*this)(i-1,j-1,k-1)
// 		   +(*this)(i+1,j+1,k-1)+(*this)(i-1,j+1,k-1)));
  return res;
}

template< >
inline __TYPE NumMatrix<__TYPE,3>::mean(const int* pos) const
{
  return mean(pos[0],pos[1],pos[2]);
}

template< >
__TYPE NumMatrix<__TYPE, 3>::
  interpol3(const int* pos, const int* diff, int r) const;

template<>
NumMatrix<__TYPE,3>& NumMatrix<__TYPE,3>::operator+=(const NumBoundary<__TYPE,3>& b);

template<>
NumMatrix<__TYPE,3>& NumMatrix<__TYPE,3>::operator-=(const NumBoundary<__TYPE,3>& b);
