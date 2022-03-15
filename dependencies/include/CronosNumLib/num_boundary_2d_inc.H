#define Boundary2d_forall(b, a) \
  {int x,y,w,l0,l1,h0,h1; \
   w = b.getWidth(); \
   l0 = b.getLow(0); \
   l1 = b.getLow(1); \
   h0 = b.getHigh(0); \
   h1 = b.getHigh(1); \
   \
   for (y=l1+w; y <= h1-w; y++) { \
     for (x=l0; x < l0+w; x++) { \
       a \
     } \
     for (x=h0-w+1; x <= h0; x++) { \
       a \
     } \
   } \
   for (x=l0; x <= h0; x++) { \
     for (y=l1; y < l1+w; y++) { \
       a \
     } \
     for (y=h1-w+1; y <= h1; y++) { \
       a \
     } \
   }} \

// ------------------------------------------------------------
// inline functions

template<>
inline __TYPE& NumBoundary<__TYPE,2>::operator()(int x, int y)
{
  assert(x >= lo[0]);
  assert(x <= hi[0]);
  assert(y >= lo[1]);
  assert(y <= hi[1]);
  assert((x < lo[0] + width) ||
	 (y < lo[1] + width) ||
	 (x > hi[0] - width) ||
	 (y > hi[1] - width));
  
  return ((__TYPE **)index)[(y-lo[1]+width)%dims[1]]
                         [(x-lo[0]+width)%dims[0]];
}

template<>
inline __TYPE NumBoundary<__TYPE,2>::operator()(int x, int y) const
{

  assert(x >= lo[0]);
  assert(x <= hi[0]);
  assert(y >= lo[1]);
  assert(y <= hi[1]);
  assert((x < lo[0] + width) ||
	 (y < lo[1] + width) ||
	 (x > hi[0] - width) ||
	 (y > hi[1] - width));
  
  return ((__TYPE **)index)[(y-lo[1]+width)%dims[1]]
                         [(x-lo[0]+width)%dims[0]];
}

template<>
inline __TYPE& NumBoundary<__TYPE,2>::operator[](const int* p)
{
  return operator()(p[0],p[1]); 
}

template<>
inline __TYPE NumBoundary<__TYPE,2>::operator[](const int* p) const
{
  return operator()(p[0],p[1]); 
}

template< >
void NumBoundary<__TYPE,2>::clear();

template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator=(const NumBoundary<__TYPE,2>& m);

template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator+=(const NumMatrix<__TYPE,2>& m);

template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator-=(const NumMatrix<__TYPE,2>& m);
				 
template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator*=(__TYPE factor);

template< >
NumBoundary<__TYPE,2>& NumBoundary<__TYPE,2>::operator/=(__TYPE factor);
				 
template< >
void NumBoundary<__TYPE,2>::deleteData();

template< >
void NumBoundary<__TYPE,2>::newData(const int* l, const int* h, int w);

template< >
void NumBoundary<__TYPE,2>::clear();

				 

