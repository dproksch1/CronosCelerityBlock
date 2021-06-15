// specific 2D methods for matrix
// -----------------------------------------------------------------------

#ifdef __sgi
template<>
ostream& ::operator<< (ostream& os, const Matrix<__TYPE, 2>& matr) 
#else
template<>
ostream& operator<< (ostream& os, const Matrix<__TYPE, 2>& matr) 
#endif
{
  if (!matr) {
    os << "(undef)" << endl;
    return os;
  }
  Index p(2);
  for (p[1] = matr.getLow()[1]; p[1] <= matr.getHigh()[1]; p[1]++) {
    for (p[0] = matr.getLow()[0]; p[0] <= matr.getHigh()[0]; p[0]++) {
      os << matr[p] << '\t';
    }
    os << endl;
  }
  return os;
}

template<>
void Matrix<__TYPE, 2>::deleteData()
{
  if (matr)
    delete[] matr;
  matr = 0;
  size = 0;

  __TYPE **_index = (__TYPE **) index;
  if (_index)
    delete[] (_index+lo[1]);
  index=0;
}

template<>
void Matrix<__TYPE, 2>::newData(const int* l, const int* h)
{
  size = 1;
  int d;
  for (d = 0; d < 2; d++) {
    assert(l[d] <= h[d]);
    lo[d] = l[d];
    hi[d] = h[d]; 
    dims[d] = h[d] - l[d] + 1;
    size *= dims[d];
  }
  matr = new __TYPE[size];

  matr_fast = matr -lo[1]*dims[0] -lo[0];

  __TYPE **_index =new __TYPE*[dims[1]];
  
  _index -= lo[1];
  _index[lo[1]] = matr - lo[0];
  
  for(int i=lo[1]+1; i<lo[1]+dims[1]; i++)
    _index[i] = _index[i-1] + dims[0];

  index=(void *)_index;
}
