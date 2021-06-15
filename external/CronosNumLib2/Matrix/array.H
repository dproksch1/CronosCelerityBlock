#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <assert.h>
#include <vector>


// ---------------------------------------------------------
// declarations

template<class T>
class Array;

/** Array of fixed size */
template<class T>
class Array {
protected:
  T* arr;
  int len;

public:
  /** constructs an Array of length len and copies the
      contents of arr into it if provided
      @memo constructor */
  Array(int len, const T* arr = 0);

  /** copy constructor */
  Array(const Array<T>&);

	/** constructs an array of zero length
	 */
	Array();

	/** resize the array to lenth len
	 */
	void resize(int len);

  /** constructs Array of length 1 containing {t0}
      @memo named constructor */
  static Array set (T t0);

  /** constructs Array of length 2 containing {t0,t1}
      @memo named constructor */
  static Array set (T t0, T t1);

  /** constructs Array of length 3 containing {t0,t1,t2}
      @memo named constructor */
  static Array set (T t0, T t1, T t2);

  /** constructs Array of length 4 containing {t0,t1,t2,t3}
      @memo named constructor */
  static Array set (T t0, T t1, T t2, T t3);

  /** constructs Array of length 5 containing {t0,t1,t2,t3,t4}
      @memo named constructor */
  static Array set (T t0, T t1, T t2, T t3, T t4);

  /** destructor */
  ~Array();

  /** returns lenght of the Array */
  int getLength() const;

  /** returns pointer to the Array, modifying forbidden */
  const T* data() const;

  /** returns pointer to the Array, modifying allowed */
  T* data();

  /** cast into C Array, modifying forbidden */
  operator const T*() const;

  /** cast into C Array, modifying allowed */
  operator T*();

  /** index operator for writing
      @memo writing */
  T& operator[](int pos);
	
	/** index operator for writing
	    @memo writing */
	T& operator()(int pos);

  /** index operator for reading
      @memo reading */
  T operator[](int pos) const;

	/** index operator for reading
	    @memo reading */
	T operator()(int pos) const;

  /** assign from std::vector */
  Array<T>& operator=(const std::vector<T>&);
  
  /** assign */
  Array<T>& operator=(const Array<T>&);

  /** comparison */
  int operator==(const Array<T>&) const;

  /** iterators */
  T* begin() { return len>0 ? &arr[0] : nullptr; }
  T* end()   { return len>0 ? &arr[len] : nullptr; }
  const T* begin() const { return len>0 ? &arr[0] : nullptr; }
  const T* end()   const { return len>0 ? &arr[len] : nullptr; }

};

/**  puts reasonable output for Array in ostream

  Example:
  \begin{verbatim}
  Array a.set(2,5);
  cout << a;
  \end{verbatim}
  @memo output for array */

template<class T>
std::ostream& operator<< (std::ostream&, const Array<T>&);

// --------------------------------------------

/** Array of fixed size which can do elementary mathematical operations */
template<class T>
class NumArray : public Array<T> {
public:
  /** constructs an NumArray of length len and copies the
      contents of arr into it if provided
      @memo constructor */
  NumArray(int len, const T* arr = 0);

  /** copy constructor */
  NumArray(const Array<T>&);


	/** empty constructor */
	NumArray();

  /** constructs NumArray of length 1 containing {t0}
      @memo named constructor */
  static NumArray set (T);

  /** constructs NumArray of length 2 containing {t0,t1}
      @memo named constructor */
  static NumArray set (T, T);

  /** constructs NumArray of length 3 containing {t0,t1,t2}
      @memo named constructor */
  static NumArray set (T, T, T);

  /** constructs NumArray of length 4 containing {t0,t1,t2,t3}
      @memo named constructor */
  static NumArray set (T t0, T t1, T t2, T t3);

  /** constructs NumArray of length 5 containing {t0,t1,t2,t3,t4}
      @memo named constructor */
  static NumArray set (T t0, T t1, T t2, T t3, T t4);

  /** constructs a linear grid from xb to xe with nx points
      @memo named constructor */
  static NumArray set_linear (T xb, T xe, int nx);
  
  /** constructs a logarithmic grid from xb to xe with nx points
    @memo named constructor */
  static NumArray set_log (T xb, T xe, int nx);

  /** set every element to zero
      @memo clear the NumArray */
  void clear();

  /** */
  NumArray<T>& operator += (const NumArray<T>&);
  /** */
  NumArray<T>& operator -= (const NumArray<T>&);
  /** */
  NumArray<T>& operator += (T);
  /** */
  NumArray<T>& operator -= (T);
  /** */
  NumArray<T>& operator *= (T);
  /** */
  NumArray<T>& operator /= (T);
  /** */
  NumArray<T> operator+ (const NumArray<T>&) const;
  /** */
  NumArray<T> operator- (const NumArray<T>&) const;
  /** */
  NumArray<T> operator+ (T) const;
  /** */
  NumArray<T> operator- (T) const;
  /** */
  NumArray<T> operator* (T) const;
  /** */
  NumArray<T> operator/ (T) const;

  /** assign from std::vector */
  NumArray<T>& operator=(const std::vector<T>&);

};

// -------------------------------------------------------------
// inline functions

template<class T>
inline T& Array<T>::operator[](int pos)
{
  return arr[pos];
}

template<class T>
inline T& Array<T>::operator()(int pos)
{
	assert(pos >= 0);
	assert(pos < len);
	return arr[pos];
}


template<class T>
inline T Array<T>::operator[](int pos) const
{
  return arr[pos];
}

template<class T>
inline T Array<T>::operator()(int pos) const
{
	assert(pos >= 0);
	assert(pos < len);
	return arr[pos];
}

template<class T>
inline const T* Array<T>::data() const
{
  return arr;
}

template<class T>
inline T* Array<T>::data()
{
  return arr;
}

template<class T>
inline Array<T>::operator const T*() const
{
  return data();
}

template<class T>
inline Array<T>::operator T*()
{
  return data();
}

/** typedef NumArray<int> Index
    @memo Integer NumArray */
typedef NumArray<int> Index;

#endif
