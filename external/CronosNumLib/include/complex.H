#ifndef COMPLEX_H
#define COMPLEX_H

#include <complex>

#if ! (__GNUC__ == 2)
typedef std::complex<double> double_complex;
#endif

#endif
