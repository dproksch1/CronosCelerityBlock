//#include "fftmatrix.H"
//
//template<class T, int rank>
//FFTMatrix<T, rank>::FFTMatrix() : 
//  NumMatrix<T, rank>(), cr_plan(0), rc_plan(0) {}
//
//template<class T, int rank>
//FFTMatrix<T, rank>::FFTMatrix(const int* d) : 
//  NumMatrix<T, rank>(d), cr_plan(0), rc_plan(0) {}
//
//template<class T, int rank>
//FFTMatrix<T, rank>::FFTMatrix(const int* l, const int* h)
//  : NumMatrix<T, rank>(l, h), cr_plan(0), rc_plan(0) {}
//
//template<class T, int rank>
//FFTMatrix<T, rank>::FFTMatrix(const Matrix<T, rank>& _matr)
//  : NumMatrix<T, rank>(_matr), cr_plan(0), rc_plan(0) {}
//
//template<class T, int rank>
//void FFTMatrix<T, rank>::newData(const int* l, const int* h)
//{
//  NumMatrix<T, rank>::newData(l, h);
//  assert(dims[0] % 2 == 0);
//  rc_plan = 
//    rfftw2d_create_plan(dims[0]-2, dims[1], FFTW_FORWARD,
//			FFTW_ESTIMATE|FFTW_IN_PLACE);
//  cr_plan = 
//    rfftw2d_create_plan(dims[0]-2, dims[1], FFTW_BACKWARD,
//			FFTW_ESTIMATE|FFTW_IN_PLACE);
//}
//
//template<class T, int rank>
//void FFTMatrix<T, rank>::deleteData()
//{
//  NumMatrix<T, rank>::deleteData();
//  if (rc_plan){
//    rfftwnd_destroy_plan(rc_plan);
//    rc_plan = 0;
//  }
//  if (cr_plan) {
//    rfftwnd_destroy_plan(cr_plan);
//    cr_plan = 0;
//  }
//}
//
//void FFTMatrix<double, 2>::crFFT()
//{
//  rfftwnd_one_complex_to_real(cr_plan, (FFTW_COMPLEX*) data(), 0);
//}
//
//void FFTMatrix<double, 2>::rcFFT()
//{
//  rfftwnd_one_real_to_complex(rc_plan, data(), 0);
//}
//
//#define instFFTMatrix(T, rank) \
//template class FFTMatrix<T, rank>; \
//
//instFFTMatrix(double, 2);
//
//#undef instFFTMatrix
//
//
