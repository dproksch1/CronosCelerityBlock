#define Boundary1d_forall(b, a) \
  {int x,w,l0,h0; \
   w = b.getWidth(); \
   l0 = b.getLow(0); \
   h0 = b.getHigh(0); \
   \
   for (x=l0; x < l0+w; x++) { \
     a \
   } \
   for (x=h0-w+1; x <= h0; x++) { \
     a \
   }} \

// ------------------------------------------------------------
// inline functions

template<>
inline __TYPE& NumBoundary<__TYPE,1>::operator()(int x)
{
  assert(x >= lo[0]);
  assert(x <= hi[0]);
  assert((x < lo[0] + width) || (x > hi[0] - width));
  
  return ((__TYPE *)index)[(x-lo[0]+width)%dims[0]];
}

template<>
inline __TYPE NumBoundary<__TYPE,1>::operator()(int x) const
{
  assert(x >= lo[0]);
  assert(x <= hi[0]);
  assert((x < lo[0] + width) || (x > hi[0] - width));
  
  return ((__TYPE *)index)[(x-lo[0]+width)%dims[0]];
}

template<>
inline __TYPE& NumBoundary<__TYPE,1>::operator[](const int* p)
{
  return operator()(p[0]); 
}

template<>
inline __TYPE NumBoundary<__TYPE,1>::operator[](const int* p) const
{
  return operator()(p[0]); 
}

template< >
void NumBoundary<__TYPE,1>::deleteData();

template< >
void NumBoundary<__TYPE,1>::newData(const int* l, const int* h, int w);

template< >
void NumBoundary<__TYPE,1>::clear();
				 
template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator=(const NumBoundary<__TYPE,1>& m);
				 
template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator+=(const NumMatrix<__TYPE,1>& m);
				 
template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator-=(const NumMatrix<__TYPE,1>& m);
				 
template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator*=(__TYPE factor);

template< >
NumBoundary<__TYPE,1>& NumBoundary<__TYPE,1>::operator/=(__TYPE factor);

