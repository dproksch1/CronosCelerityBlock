#ifndef CRONOS_VECTOR_H
#define CRONOS_VECTOR_H = 1

namespace cronos {
	template <class T>
	class vector {
	public:
		vector(){
			val[0] = 0;
			val[1] = 0;
			val[2] = 0;
		};
		vector(T x, T y, T z){
			this->val[0] = x;
			this->val[1] = y;
			this->val[2] = z;
		};
		T abs();
		T quad();
		void populate(T x, T y, T z);
		T get(const int &);
		void set(const int &, const T &);
	private:
		T val[3];
	};

	template <class T> 
	void vector<T>::populate(T x, T y, T z)
	{
		val[0] = x;
		val[1] = y;
		val[2] = z;
	}
  
	template <class T>
	T vector<T>::get(const int &comp) 
	{
		return val[comp];
	}

	template <class T>
	void vector<T>::set(const int &comp, const T &val)
	{
		this->val[comp] = val;
	}
  
	template <class T> 
	T vector<T>::abs()
	{
		return sqrt(sqr(val[0]) + sqr(val[1]) + sqr(val[2]));
	}
  
	template <class T> 
	T vector<T>::quad()
	{
		return sqr(val[0]) + sqr(val[1]) + sqr(val[2]);
	}
}//namespace cronos

#endif
