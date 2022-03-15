#ifndef CRONOS_TEST_DATA_H
#define CRONOS_TEST_DATA_H = 1

#include <matrix.H>

//#define REAL double

#ifndef DATA_TYPE
#define DATA_TYPE double
#define REAL_FLAG
#endif

#ifndef DIM
#define DIM 3
#endif

class Pot : public NumMatrix<DATA_TYPE,DIM>
{
public:
    Pot();
};

#endif