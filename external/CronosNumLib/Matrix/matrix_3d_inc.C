// specific 3D methods for matrix
// -----------------------------------------------------------------------

#ifdef __sgi
template<>
ostream& ::operator<< (ostream& os, const Matrix<__TYPE, 3>& matr) 
#else
template<>
ostream& operator<< (ostream& os, const Matrix<__TYPE, 3>& matr) 
#endif
{
  if (!matr) {
    os << "(undef)" << endl;
    return os;
  }
  Index p(3);
  for (p[2] = matr.getLow()[2]; p[2] <= matr.getHigh()[2]; p[2]++) {
    for (p[1] = matr.getLow()[1]; p[1] <= matr.getHigh()[1]; p[1]++) {
      for (p[0] = matr.getLow()[0]; p[0] <= matr.getHigh()[0]; p[0]++) {
	os << matr[p] << '\t';
      }
      os << endl;
    }
    os << endl;
  }
  return os;
}

template<>
void Matrix<__TYPE, 3>::deleteData()
{
  if (matr)
    delete[] matr;
  matr = 0;
  size = 0;

  __TYPE ***_index = (__TYPE ***) index;
  if (_index) {
    for (int z=lo[2]; z < lo[2]+dims[2]; z++)
      delete[] (_index[z]+lo[1]);
    delete[] (_index+lo[2]);
  }

  index=0;
}

template<>
void Matrix<__TYPE, 3>::newData(const int* l, const int* h)
{
  size = 1;
  int d;
  for (d = 0; d < 3; d++) {
    assert(l[d] <= h[d]);
    lo[d] = l[d];
    hi[d] = h[d]; 
    dims[d] = h[d] - l[d] + 1;
    size *= dims[d];
  }
  matr = new __TYPE[size];

  matr_fast = matr -((lo[2]*dims[1]+lo[1])*dims[0]+lo[0]);

  __TYPE ***_index =new __TYPE**[dims[2]];
  
  _index -= lo[2];
  for (int z=lo[2]; z < lo[2]+dims[2]; z++) {
    _index[z] = new __TYPE*[dims[1]];
    _index[z] -= lo[1];
    _index[z][lo[1]] = matr+((z-lo[2])*dims[0]*dims[1]) - lo[0];

    for (int y=lo[1]+1; y<lo[1]+dims[1]; y++)
      _index[z][y] = _index[z][y-1] + dims[0];
  }

  index=(void *)_index;
}
