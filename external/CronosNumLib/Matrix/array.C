#include "array.H"
#include <cmath>

using namespace std;

// --------------------------------------------------------------
// implementation

template<class T>
Array<T>::Array(int _len, const T* _arr)
{
  arr = new T[len = _len];
  if (_arr)
    for (int i = 0; i < len; i++)
      arr[i] = _arr[i];
}

template<class T>
Array<T>::Array(const Array<T>& _arr)
{
  arr = new T[len = _arr.len];
  for (int i = 0; i < len; i++)
    arr[i] = _arr[i];
}

template<class T>
Array<T>::Array() {
	len = 0;
	arr = NULL;
}

template<class T>
void Array<T>::resize(int _len) {
	if(arr != NULL) {
		delete[] arr;
	}
	arr = new T[len = _len];
	for (int i = 0; i < this->len; i++) {
		this->arr[i] = 0;
	}
}

template<class T>
Array<T> Array<T>::set(T t0)
{
  T tmp[1] = { t0 };
  return Array(1, tmp);
}

template<class T>
Array<T> Array<T>::set(T t0, T t1)
{
  T tmp[2] = { t0, t1 };
  return Array(2, tmp);
}

template<class T>
Array<T> Array<T>::set(T t0, T t1, T t2)
{
  T tmp[3] = { t0, t1, t2 };
  return Array(3, tmp);
}

template<class T>
Array<T> Array<T>::set(T t0, T t1, T t2, T t3)
{
  T tmp[4] = { t0, t1, t2, t3 };
  return Array(4, tmp);
}

template<class T>
Array<T> Array<T>::set(T t0, T t1, T t2, T t3, T t4)
{
  T tmp[5] = { t0, t1, t2, t3, t4 };
  return Array(5, tmp);
}


template<class T>
Array<T>::~Array()
{
	if(arr != NULL) {
		delete[] arr;
	}
}

template<class T>
int Array<T>::getLength() const
{
  return len;
}

// template<class T>
// Array<T>& Array<T>::operator=(const Array<T>& _arr)
// {
//   assert(len == _arr.len);
//   for (int i = 0; i < len; i++)
//     arr[i] = _arr[i];
//   return *this;
// }

template<class T>
Array<T>& Array<T>::operator=(const Array<T>& _arr)
{
	if(len != _arr.len) {
		arr = new T[len = _arr.len];
	}
	for (int i = 0; i < len; i++)
		arr[i] = _arr[i];
	return *this;
}

template<class T>
Array<T>& Array<T>::operator=(const std::vector<T>& _vec)
{
  if(len != _vec.size()) {
    arr = new T[len = _vec.size()];
  }
  for (int i = 0; i < len; i++)
    arr[i] = _vec[i];
  return *this;
}

template<class T>
int Array<T>::operator==(const Array<T>& _arr) const
{
  if (len != _arr.len)
    return 0;

  for (int i = 0; i < len; i++)
    if (arr[i] != _arr[i])
      return 0;
  return 1;
}

template<class T>
ostream& operator<< (ostream& os, const Array<T>& arr)
{
  os << "[ ";
  for (int i = 0; i < arr.getLength(); i++) {
    if (i > 0)
      os << " ";
    os << arr[i];
  }
  os << " (" << arr.getLength() << ")]";
  return os;
}

// --------------------------------------

template<class T>
NumArray<T>::NumArray(int _len, const T* _arr)
  : Array<T>(_len, _arr) {}

template<class T>
NumArray<T>::NumArray(const Array<T>& _arr) : Array<T>(_arr) {}

template<class T>
NumArray<T>::NumArray() : Array<T>() {}

template<class T>
NumArray<T> NumArray<T>::set(T t0)
{
  T tmp[1] = { t0 };
  return NumArray(1, tmp);
}

template<class T>
NumArray<T> NumArray<T>::set(T t0, T t1)
{
  T tmp[2] = { t0, t1 };
  return NumArray(2, tmp);
}

template<class T>
NumArray<T> NumArray<T>::set(T t0, T t1, T t2)
{
  T tmp[3] = { t0, t1, t2 };
  return NumArray(3, tmp);
}

template<class T>
NumArray<T> NumArray<T>::set(T t0, T t1, T t2, T t3)
{
  T tmp[4] = { t0, t1, t2, t3 };
  return NumArray(4, tmp);
}

template<class T>
NumArray<T> NumArray<T>::set(T t0, T t1, T t2, T t3, T t4)
{
  T tmp[5] = { t0, t1, t2, t3, t4 };
  return NumArray(5, tmp);
}

template<class T>
NumArray<T> NumArray<T>::set_linear (T xb, T xe, int nx)
{
	NumArray tmp(nx);
	tmp(0) = xb;
	tmp(nx-1) = xe;
	T dx = (xe - xb) / (nx - 1);
	for(int ix = 1; ix < nx - 1; ++ix) {
		tmp(ix) = tmp(ix - 1) + dx;
	}
  return tmp;
}

template<class T>
NumArray<T> NumArray<T>::set_log (T xb, T xe, int nx)
{
	NumArray tmp(nx);
	tmp(0) = xb;
	tmp(nx-1) = xe;
	T dx = pow(xe / xb, 1. / (nx - 1));
	for(int ix = 1; ix < nx - 1; ++ix) {
		tmp(ix) = tmp(ix - 1) * dx;
	}
  return tmp;
}

template<class T>
void NumArray<T>::clear()
{
  for (int i = 0; i < this->len; i++)
    this->arr[i] = 0;
}

template<class T>
NumArray<T>& NumArray<T>::operator*=(T t)
{
  for (int i = 0; i < this->len; i++)
    this->arr[i] *= t;
  return *this;
}

template<class T>
NumArray<T>& NumArray<T>::operator/=(T t)
{
  for (int i = 0; i < this->len; i++)
    this->arr[i] /= t;
  return *this;
}

template<class T>
NumArray<T>& NumArray<T>::operator+=(const NumArray<T>& _arr)
{
  assert(this->len == _arr.len);
  for (int i = 0; i < this->len; i++)
    this->arr[i] += _arr[i];
  return *this;
}

template<class T>
NumArray<T>& NumArray<T>::operator-=(const NumArray<T>& _arr)
{
  assert(this->len == _arr.len);
  for (int i = 0; i < this->len; i++)
    this->arr[i] -= _arr[i];
  return *this;
}

template<class T>
NumArray<T>& NumArray<T>::operator+=(T t)
{
  for (int i = 0; i < this->len; i++)
    this->arr[i] += t;
  return *this;
}

template<class T>
NumArray<T>& NumArray<T>::operator-=(T t)
{
  for (int i = 0; i < this->len; i++)
    this->arr[i] -= t;
  return *this;
}

template<class T>
NumArray<T> NumArray<T>::operator+(const NumArray<T>& _arr) const
{
  assert(this->len == _arr.len);
  NumArray<T> tmp = *this;
  tmp += _arr;
  return tmp;
}

template<class T>
NumArray<T> NumArray<T>::operator-(const NumArray<T>& _arr) const
{
  assert(this->len == _arr.len);
  NumArray<T> tmp = *this;
  tmp -= _arr;
  return tmp;
}

template<class T>
NumArray<T> NumArray<T>::operator+(T t) const
{
  NumArray<T> tmp = *this;
  tmp += t;
  return tmp;
}

template<class T>
NumArray<T> NumArray<T>::operator-(T t) const
{
  NumArray<T> tmp = *this;
  tmp -= t;
  return tmp;
}

template<class T>
NumArray<T> NumArray<T>::operator*(T t) const
{
  NumArray<T> tmp = *this;
  tmp *= t;
  return tmp;
}

template<class T>
NumArray<T> NumArray<T>::operator/(T t) const
{
  NumArray<T> tmp = *this;
  tmp /= t;
  return tmp;
}


template<class T>
NumArray<T>& NumArray<T>::operator=(const std::vector<T>& _vec)
{
  if(this->len != _vec.size()) {
    this->arr = new T[this->len = _vec.size()];
  }
  for (int i = 0; i < this->len; i++)
    this->arr[i] = _vec[i];
  return *this;
}

// ------------------------------------------------------------
// explicit instantiation

#ifndef __sgi
#define instNumArray(T) \
template class Array<T>; \
template class NumArray<T>; \
template ostream& operator<<(ostream&, const Array<T> &);
#else
template class Array<T>; \
template class NumArray<T>; \
template ostream& ::operator<<(ostream&, const Array<T> &);
#endif

instNumArray(int)
instNumArray(long)
instNumArray(float)
instNumArray(double)

#undef instNumArray


