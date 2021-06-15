template< >
void NumBoundary<__TYPE,3>::deleteData()
{
  __TYPE ***index_ = (__TYPE ***)index;
  if (index_) {
    for (int z=0; z < dims[2]; z++) {
      for(int y=0; y < dims[1]; y++)
	delete[] index_[z][y];
      delete[] index_[z];
    }
    delete[] index_;
  }
  index = 0;
}

template< >
void NumBoundary<__TYPE,3>::newData(const int* l, const int* h, int w)
{
  width = w;
  for (int d = 0; d < 3; d++) {
    assert(l[d] <= h[d]);
    lo[d] = l[d];
    hi[d] = h[d]; 
    dims[d] = h[d] - l[d] + 1;
  }
 
  int y,z;

  __TYPE ***index_ = new __TYPE**[dims[2]];

  for(z = 0; z < 2*width; z++) {
    index_[z] = new __TYPE*[dims[1]];
    for (y = 0; y < dims[1]; y++)
      index_[z][y] = new __TYPE[dims[0]];
  }

  for(; z < dims[2]; z++) {
    index_[z] = new __TYPE*[dims[1]];
    for(y=0; y < 2*width; y++)
      index_[z][y] = new __TYPE[dims[0]];
    
    for(; y < dims[1]; y++)
      index_[z][y] = new __TYPE[2*width];
  }

  index = (void *) index_;
}

template< >
void NumBoundary<__TYPE,3>::clear()
{
  Boundary3d_forall((*this), (*this)(x,y,z) = 0;)
}
				 
template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator=(const NumBoundary<__TYPE,3>& m)
{
  Boundary3d_forall((*this), (*this)(x,y,z) = m(x,y,z);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator+=(const NumMatrix<__TYPE,3>& m)
{
  Boundary3d_forall((*this), (*this)(x,y,z) += m(x,y,z);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator-=(const NumMatrix<__TYPE,3>& m)
{
  Boundary3d_forall((*this), (*this)(x,y,z) -= m(x,y,z);)
  return (*this);
}
				 
template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator*=(__TYPE factor)
{
  Boundary3d_forall((*this), (*this)(x,y,z) *= factor;)
  return(*this);  
}

template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator/=(__TYPE factor)
{
  __TYPE inv = __TYPE(1./factor);
  Boundary3d_forall((*this), (*this)(x,y,z) *= inv;)
  return(*this);  
}

