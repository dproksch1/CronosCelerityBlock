#include "matrix_2d_inc.H"

template< >
inline __TYPE NumMatrix<__TYPE,2>::mean(int i, int j) const
{
  __TYPE res = __TYPE(0.0625*((*this)(i+1,j+1) +(*this)(i+1,j-1)
		  +(*this)(i-1,j+1)+(*this)(i-1,j-1))
    +0.125*((*this)(i,j+1)+(*this)(i,j-1)
	    +(*this)(i+1,j)+(*this)(i-1,j)) 
    +0.25*(*this)(i,j));
  return res;
}


// // #if (__TYPE == double_complex)
// template< >
// inline __TYPE NumMatrix<__TYPE,2>::get_max() const
// {
// 	__TYPE maximum = this->matr[0];
// 	for (int i = 0; i < this->size; i++) {
// 		maximum = std::max(maximum, this->matr[i]);
// 	}
// 	return maximum;
// }

// template< >
// inline __TYPE NumMatrix<__TYPE,2>::get_min() const
// {
// 	__TYPE minimum = this->matr[0];
// 	for (int i = 0; i < this->size; i++) {
// 		minimum = std::min(minimum, this->matr[i]);
// 	}
// 	return minimum;
// }

// // #endif


template< >
inline __TYPE NumMatrix<__TYPE,2>::mean(const int* pos) const
{
  return mean(pos[0],pos[1]);
}

template< >
__TYPE NumMatrix<__TYPE, 2>::
  interpol3(const int* pos, const int* diff, int r) const;

template<>
NumMatrix<__TYPE,2>& NumMatrix<__TYPE,2>::operator+=(const NumBoundary<__TYPE,2>& b);

template<>
NumMatrix<__TYPE,2>& NumMatrix<__TYPE,2>::operator-=(const NumBoundary<__TYPE,2>& b);


