template< >
void NumBoundary<__TYPE,2>::deleteData()
{
  __TYPE **index_ = (__TYPE **)index;
  if (index_) {
    for(int y=0; y < dims[1]; y++)
      delete[] index_[y];
    delete[] index_;
  }
  index = 0;
}

template< >
void NumBoundary<__TYPE,2>::newData(const int* l, const int* h, int w)
{
  width = w;
  for (int d = 0; d < 2; d++) {
    assert(l[d] <= h[d]);
    lo[d] = l[d];
    hi[d] = h[d]; 
    dims[d] = h[d] - l[d] + 1;
  }
 
  int y;

  __TYPE **index_ = new __TYPE*[dims[1]];
  
  for(y = 0; y < 2*width; y++)
    index_[y] = new __TYPE[dims[0]];
  
  for(; y < dims[1]; y++)
    index_[y] = new __TYPE[2*width];

  index = (void *) index_;
}

template< >
void NumBoundary<__TYPE,2>::clear()
{
  Boundary2d_forall((*this), (*this)(x,y) = 0;)
}
				 
template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator=(const NumBoundary<__TYPE,2>& m)
{
  Boundary2d_forall((*this), (*this)(x,y) = m(x,y);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator+=(const NumMatrix<__TYPE,2>& m)
{
  Boundary2d_forall((*this), (*this)(x,y) += m(x,y);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator-=(const NumMatrix<__TYPE,2>& m)
{
  Boundary2d_forall((*this), (*this)(x,y) -= m(x,y);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator*=(__TYPE factor)
{
  Boundary2d_forall((*this), (*this)(x,y) *= factor;)
  return(*this);  
}

template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator/=(__TYPE factor)
{
  __TYPE inv = __TYPE(1./factor);
  Boundary2d_forall((*this), (*this)(x,y) *= inv;)
  return(*this);  
}


