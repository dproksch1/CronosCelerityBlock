#define Boundary3d_forall(b, a) \
  {int x,y,z,w,l0,l1,l2,h0,h1,h2; \
   w = b.getWidth(); \
   l0 = b.getLow(0); \
   l1 = b.getLow(1); \
   l2 = b.getLow(2); \
   h0 = b.getHigh(0); \
   h1 = b.getHigh(1); \
   h2 = b.getHigh(2); \
   \
   for (z=l2+w; z <= h2-w; z++) { \
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
     } \
   } \
   for (y=l1; y <= h1; y++) { \
     for (x=l0; x <= h0; x++) { \
       for (z=l2; z < l2+w; z++) { \
	 a \
       } \
       for (z=h2-w+1; z <= h2; z++) { \
	 a \
       } \
     } \
   }} \

// ------------------------------------------------------------
// inline functions

template<>
inline __TYPE& NumBoundary<__TYPE,3>::operator()(int x, int y, int z)
{
  assert(x >= lo[0]);
  assert(x <= hi[0]);
  assert(y >= lo[1]);
  assert(y <= hi[1]);
  assert(z >= lo[2]);
  assert(z <= hi[2]);
  assert((x < lo[0] + width) ||
	 (y < lo[1] + width) ||
	 (z < lo[2] + width) ||
	 (x > hi[0] - width) ||
	 (y > hi[1] - width) ||
	 (z > hi[2] - width));
  
  return ((__TYPE ***)index)[(z-lo[2]+width)%dims[2]]
                          [(y-lo[1]+width)%dims[1]]
                          [(x-lo[0]+width)%dims[0]];
}

template<>
inline __TYPE NumBoundary<__TYPE,3>::operator()(int x, int y, int z) const
{
  assert(x >= lo[0]);
  assert(x <= hi[0]);
  assert(y >= lo[1]);
  assert(y <= hi[1]);
  assert(z >= lo[2]);
  assert(z <= hi[2]);
  assert((x < lo[0] + width) ||
	 (y < lo[1] + width) ||
	 (z < lo[2] + width) ||
	 (x > hi[0] - width) ||
	 (y > hi[1] - width) ||
	 (z > hi[2] - width));
  
  return ((__TYPE ***)index)[(z-lo[2]+width)%dims[2]]
                          [(y-lo[1]+width)%dims[1]]
                          [(x-lo[0]+width)%dims[0]];
}

template<>
inline __TYPE& NumBoundary<__TYPE,3>::operator[](const int* p)
{
  return operator()(p[0],p[1],p[2]); 
}

template<>
inline __TYPE NumBoundary<__TYPE,3>::operator[](const int* p) const
{
  return operator()(p[0],p[1],p[2]); 
}

  
template< >
void NumBoundary<__TYPE,3>::deleteData();

template< >
void NumBoundary<__TYPE,3>::newData(const int* l, const int* h, int w);

template< >
void NumBoundary<__TYPE,3>::clear();
				 
template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator=(const NumBoundary<__TYPE,3>& m);
				 
template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator+=(const NumMatrix<__TYPE,3>& m);
				 
template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator-=(const NumMatrix<__TYPE,3>& m);
				 
template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator*=(__TYPE factor);

template< >
NumBoundary<__TYPE,3>& NumBoundary<__TYPE,3>::operator/=(__TYPE factor);

