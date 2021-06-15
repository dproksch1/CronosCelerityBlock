template< >
inline  __TYPE& Matrix<__TYPE, 2>::operator ()(int i, int j)
{
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]) ;
  assert(j <= hi[1]);
  return ((__TYPE**)index)[j][i];
}

template< >
inline  __TYPE Matrix<__TYPE, 2>::operator ()(int i, int j) const
{
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  return ((__TYPE**)index)[j][i];
}

template< >
inline  __TYPE& Matrix<__TYPE, 2>::operator[](const int* pos)
{
  return operator ()(pos[0],pos[1]);
}

template< >
inline  __TYPE Matrix<__TYPE, 2>::operator[](const int* pos) const
{
  return operator ()(pos[0],pos[1]);
}

template<>
void Matrix<__TYPE, 2>::deleteData();

template<>
void Matrix<__TYPE, 2>::newData(const int* l, const int* h);

#ifdef __sgi
template<>
std::ostream& ::operator<< (std::ostream&, const Matrix<__TYPE, 2>&);
#else
template<>
std::ostream& operator<< (std::ostream&, const Matrix<__TYPE, 2>&);
#endif
