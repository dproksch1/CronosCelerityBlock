template< >
inline __TYPE& Matrix<__TYPE, 4>::operator ()(int i, int j, int k, int l)
{
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  assert(k >= lo[2]);
  assert(k <= hi[2]);
  assert(l >= lo[3]);
  assert(l <= hi[3]);
  return matr[((((l-lo[3]))*dims[2] + 
                (k-lo[2]))*dims[1] + (j-lo[1]))*dims[0] + (i - lo[0])];
}

template< >
inline __TYPE Matrix<__TYPE, 4>::operator ()(int i, int j, int k, int l) const
{
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  assert(k >= lo[2]);
  assert(k <= hi[2]);
  assert(l >= lo[3]);
  assert(l <= hi[3]);
  return matr[((((l-lo[3]))*dims[2] + 
                (k-lo[2]))*dims[1] + (j-lo[1]))*dims[0] + (i - lo[0])];
}

