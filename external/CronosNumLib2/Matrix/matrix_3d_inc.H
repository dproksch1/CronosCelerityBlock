template< >
inline __TYPE& Matrix<__TYPE, 3>::operator ()(int i, int j, int k)
{
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  assert(k >= lo[2]);
  assert(k <= hi[2]);
  return ((__TYPE***)index)[k][j][i];
//  return matr_fast[i + dims[0]*(j + dims[1]*k)];
}

template< >
inline __TYPE Matrix<__TYPE, 3>::operator ()(int i, int j, int k) const
{
  assert(i >= lo[0]);
  assert(i <= hi[0]);
  assert(j >= lo[1]);
  assert(j <= hi[1]);
  assert(k >= lo[2]);
  assert(k <= hi[2]);
  return ((__TYPE***)index)[k][j][i];
}

template< >
inline __TYPE& Matrix<__TYPE, 3>::operator[](const int* pos)
{
  return operator ()(pos[0],pos[1],pos[2]);
}

template< >
inline __TYPE Matrix<__TYPE, 3>::operator[](const int* pos) const
{
  return operator ()(pos[0],pos[1],pos[2]);
}

template<>
void Matrix<__TYPE, 3>::deleteData();

template<>
void Matrix<__TYPE, 3>::newData(const int* l, const int* h);

template<>
std::ostream& operator<< (std::ostream&, const Matrix<__TYPE, 3>&);
