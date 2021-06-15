// specific 1D methods for matrix
// -----------------------------------------------------------------------

#ifdef __sgi
template<>
ostream& ::operator<< (ostream& os, const Matrix<__TYPE, 1>& matr) 
#else
template<>
ostream& operator<< (ostream& os, const Matrix<__TYPE, 1>& matr) 
#endif
{
  if (!matr) {
    os << "(undef)" << endl;
    return os;
  }
  Index p(1);
  for (p[0] = matr.getLow()[0]; p[0] <= matr.getHigh()[0]; p[0]++) {
    os << matr[p] << '\t';
  }
  os << endl;
  return os;
}
