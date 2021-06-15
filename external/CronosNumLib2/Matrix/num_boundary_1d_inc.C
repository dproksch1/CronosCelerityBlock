template< >
void NumBoundary<__TYPE,1>::deleteData()
{
  __TYPE *index_ = (__TYPE *)index;
  if (index_) {
    delete[] index_;
  }
  index = 0;
}

template< >
void NumBoundary<__TYPE,1>::newData(const int* l, const int* h, int w)
{
  width = w;
  assert(l[0] <= h[0]);
  lo[0] = l[0];
  hi[0] = h[0]; 
  dims[0] = h[0] - l[0] + 1;
 
  __TYPE *index_ = new __TYPE[2*width];
  
  index = (void *) index_;
}

template< >
void NumBoundary<__TYPE,1>::clear()
{
  Boundary1d_forall((*this), (*this)(x) = 0;)
}
				 
template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator=(const NumBoundary<__TYPE,1>& m)
{
  Boundary1d_forall((*this), (*this)(x) = m(x);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator+=(const NumMatrix<__TYPE,1>& m)
{
  Boundary1d_forall((*this), (*this)(x) += m(x);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator-=(const NumMatrix<__TYPE,1>& m)
{
  Boundary1d_forall((*this), (*this)(x) -= m(x);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator*=(__TYPE factor)
{
  Boundary1d_forall((*this), (*this)(x) *= factor;)
  return(*this);  
}

template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator/=(__TYPE factor)
{
  __TYPE inv = __TYPE(1./factor);
  Boundary1d_forall((*this), (*this)(x) *= inv;)
  return(*this);  
}
