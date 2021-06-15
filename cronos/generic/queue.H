#pragma once

#include <CL/sycl.hpp>
#include <matrix.H>

#include <fstream>

using Queue = cl::sycl::queue;

using DeviceSelector = cl::sycl::host_selector;

template<typename DataType, int Dims>
using Buffer = cl::sycl::buffer<DataType, Dims>;

template<int Dims>
using Range = cl::sycl::range<Dims>;

template<int Dims>
using Item = cl::sycl::item<Dims>;

template<int Dims>
using Id = cl::sycl::id<Dims>;

using Handler = cl::sycl::handler;

#define __xstr(a) __insieme_utils_str(a)
#define __insieme_utils_str(a) #a

#ifdef __GNUC__
#define __insieme_unused __attribute__((unused))
#else
#define __insieme_unused
#endif

namespace insieme {
	namespace utils {

		// allows to provide a lambda that adds printed information when an assertion occurs
		//void setAssertExtraInfoPrinter(std::function<void(void)> printer);
		// clears the extra information printer lambda
		//void clearAssertExtraInfoPrinter();

		namespace detail {

			//extern thread_local std::function<void(void)> extraAssertInformationPrinter;

			struct LazyAssertion {
				bool value;
				LazyAssertion(bool value) : value(value) {}
				~LazyAssertion() /*{ if(!value) __debugbreak(); } //*/= default;
				operator bool() const {
					return !value;
				}
			};

		} // end namespace detail
	} // end namespace utils
} // end namespace insieme

#define assert_true(_COND)                                                                                                                                     \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion((bool)(_COND)))                                                                       \
	std::cerr << "\nAssertion " #_COND " of " __FILE__ ":" __xstr(__LINE__) " failed!\n"

#define assert_eq(_A, _B)                                                                                                                                      \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion((_A) == (_B)))                                                                        \
	std::cerr << "\nAssertion " #_A " == " #_B " of " __FILE__ ":" __xstr(__LINE__) " failed!\n\t" #_A " = " << (_A) << "\n\t" #_B " = " << (_B) << "\n"

#define assert_ne(_A, _B)                                                                                                                                      \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion((_A) != (_B)))                                                                        \
	std::cerr << "\nAssertion " #_A " != " #_B " of " __FILE__ ":" __xstr(__LINE__) " failed!\n\t" #_A " = " << (_A) << "\n\t" #_B " = " << (_B) << "\n"

#define assert_lt(_A, _B)                                                                                                                                      \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion((_A) < (_B)))                                                                         \
	std::cerr << "\nAssertion " #_A " < " #_B " of " __FILE__ ":" __xstr(__LINE__) " failed!\n\t" #_A " = " << (_A) << "\n\t" #_B " = " << (_B) << "\n"

#define assert_le(_A, _B)                                                                                                                                      \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion((_A) <= (_B)))                                                                        \
	std::cerr << "\nAssertion " #_A " <= " #_B " of " __FILE__ ":" __xstr(__LINE__) " failed!\n\t" #_A " = " << (_A) << "\n\t" #_B " = " << (_B) << "\n"

#define assert_gt(_A, _B)                                                                                                                                      \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion((_A) > (_B)))                                                                         \
	std::cerr << "\nAssertion " #_A " > " #_B " of " __FILE__ ":" __xstr(__LINE__) " failed!\n\t" #_A " = " << (_A) << "\n\t" #_B " = " << (_B) << "\n"

#define assert_ge(_A, _B)                                                                                                                                      \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion((_A) >= (_B)))                                                                        \
	std::cerr << "\nAssertion " #_A " >= " #_B " of " __FILE__ ":" __xstr(__LINE__) " failed!\n\t" #_A " = " << (_A) << "\n\t" #_B " = " << (_B) << "\n"

#define assert_fail()                                                                                                                                          \
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion(false)) std::cerr << "\nAssertion failed in " __FILE__ ":" __xstr(__LINE__) " - "



enum Direction {
	DirX = 0,
	DirY = 1,
	DirZ = 2,

	DirMax = 3,
};

//template<typename DataType>
//Buffer<DataType, 1> toSYCL(const NumMatrix<DataType, 1>& matrix) {
//	const int dimIndex = 0;
//	const auto& indexStart = matrix.getLow(dimIndex);
//	const auto& indexEnd = matrix.getHigh(dimIndex);
//	const auto& size = indexEnd - indexStart + 1;
//
//	assert_ge(size, 0) << "Expected size to be greater than zero" << std::endl;
//
//	Buffer<DataType, 1> retVal = Buffer<DataType, 1>(Range<1>(size));
//	auto retValAccessor = retVal.get_access<cl::sycl::access::mode::discard_write>();
//	
//	for (int i = 0; i < size; ++i) {
//		retValAccessor[i] = matrix(i+indexStart);
//	}
//
//	return retVal;
//}

template<typename DataType>
Buffer<DataType, 1> extractPoleToSYCL(const NumMatrix<DataType, 3>& matrix, size_t size, int dir, int offsetI, int offsetJ, int offsetK) {
	Buffer<DataType, 1> retVal = Buffer<DataType, 1>(Range<1>(size));

	auto retVal_acc = retVal.get_access<cl::sycl::access::mode::discard_write>(NULL);

	//std::cout << "size: " << size << ", dir: " << dir << ", offsetI: " << offsetI << ", offsetJ: " << offsetJ << ", offsetK: " << offsetK << std::endl;

	if (dir == 0) {
		// x direction
		for (int i = 0; i < size; ++i) {
			retVal_acc[i] = matrix(i + offsetI, offsetJ, offsetK);
		}
	} else if (dir == 1) {
		// y direction
		for (int i = 0; i < size; ++i) {
			retVal_acc[i] = matrix(offsetI, i + offsetJ, offsetK);
		}
	} else if (dir == 2) {
		// z direction
		for (int i = 0; i < size; ++i) {
			retVal_acc[i] = matrix(offsetI, offsetJ, i + offsetK);
		}
	} else {
		assert_fail() << "Not implemented" << std::endl;
	}

	return retVal;
}

template<typename DataType>
NumMatrix<DataType, 1> toCronos(Buffer<DataType, 1> buffer, int indexStart, int indexEnd) {
	NumMatrix<DataType, 1> retVal;
	const auto& size = indexEnd - indexStart + 1;

	assert_ge(size, 0) << "Expected size to be greater than zero" << std::endl;
	assert_eq(size, buffer.get_count()) << "Number of elements doesn't match" << std::endl;

	retVal.resize(&indexStart, &indexEnd);

	auto bufferAccessor = buffer.get_access<cl::sycl::access::mode::read>(NULL);

	for (int i = 0; i < size; ++i) {
		retVal(i+indexStart) = bufferAccessor[i];
	}

	return retVal;
}

template<typename DataType>
bool isEqual(const NumMatrix<DataType, 1>& matrix, Buffer<DataType, 1> buffer) {
	const int dimIndex = 0;
	const auto& syclData = toCronos(buffer, matrix.getLow(dimIndex), matrix.getHigh(dimIndex));
	const auto& origData = matrix;

	// check dimension lengths
	assert_eq(origData.getHigh(dimIndex), syclData.getHigh(dimIndex));
	assert_eq(origData.getLow(dimIndex), syclData.getLow(dimIndex));

	bool equal = true;
	for (int i = syclData.getLow(dimIndex); i <= syclData.getHigh(dimIndex); ++i) {
		assert_eq(origData(i), syclData(i)) << "Mismatching elements for i = " << i << std::endl;
		if (origData(i) != syclData(i)) {
			equal = false;
			break;
		}
	}
	return equal;
}

template<typename DataType>
bool isEqual(Buffer<DataType, 1> buffer, const NumMatrix<DataType, 1>& matrix) {
	return isEqual(matrix, buffer);
}

#define assert_sycl_eq(_A, _B)                                                                                                                                      \
	queue.wait();																																					\
	if(__insieme_unused auto __x = insieme::utils::detail::LazyAssertion(isEqual((_A), (_B))))                                                                      \
	std::cerr << "\nAssertion " #_A " == " #_B " of " __FILE__ ":" __xstr(__LINE__) " failed!\n"


const std::string debugFilename = "debug.new.txt";

template<typename DataType>
void dump(const NumMatrix<DataType, 1>& matrix) {
	std::ofstream outfile;
	outfile.open(debugFilename.c_str(), std::ios_base::app);
	outfile << "dims: ";
	for (int dim = 0; dim < 1; ++dim) {
		outfile << dim << ": [" << matrix.getLow(dim) << ", " << matrix.getHigh(dim) << "], ";
	}
	outfile << std::endl << std::flush;

	for (int i = matrix.getLow(0); i <= matrix.getHigh(0); ++i) {
		outfile << "(" << i << "):" << matrix(i) << " ";
	}
	outfile << std::flush;
}

template<typename DataType>
void dump(const NumMatrix<DataType, 2>& matrix) {
	std::ofstream outfile;
	outfile.open(debugFilename.c_str(), std::ios_base::app);
	outfile << "dims: ";
	for (int dim = 0; dim < 2; ++dim) {
		outfile << dim << ": [" << matrix.getLow(dim) << ", " << matrix.getHigh(dim) << "], ";
	}
	outfile << std::endl << std::flush;
	
	for (int i = matrix.getLow(0); i <= matrix.getHigh(0); ++i) {
		for (int j = matrix.getLow(1); j <= matrix.getHigh(1); ++j) {
			outfile << "(" << i << "," << j << "):" << matrix(i,j) << " ";
		}
	}
	outfile << std::flush;
}

template<typename DataType>
void dump(const NumMatrix<DataType, 3>& matrix) {
	std::ofstream outfile;
	outfile.open(debugFilename.c_str(), std::ios_base::app);
	outfile << "dims: ";
	for (int dim = 0; dim < 3; ++dim) {
		outfile << dim << ": [" << matrix.getLow(dim) << ", " << matrix.getHigh(dim) << "], ";
	}
	outfile << std::endl << std::flush;

	for (int i = matrix.getLow(0); i <= matrix.getHigh(0); ++i) {
		for (int j = matrix.getLow(1); j <= matrix.getHigh(1); ++j) {
			for (int k = matrix.getLow(1); k <= matrix.getHigh(2); ++k) {
				outfile << "(" << i << "," << j << "," << k << "):" << matrix(i, j, k) << " ";
			}
		}
	}
	outfile << std::flush;
}

template<typename DataType>
void clearBuffer(Queue queue, Buffer<DataType, 1> buffer) {
#ifdef USE_SYCL
	queue.submit([&](sycl::handler& cgh) {
		auto acc = buffer.get_access<cl::sycl::access::mode::write>(cgh);
		size_t size = buffer.get_count();
		cgh.parallel_for<class ClearBuffer>(Range<1>(size), [=](Item<1> item) {
			size_t i = item.get(0);
			acc[i] = DataType{};
		});
	});
#endif
}

template<typename DataType>
void clearBufferOnHost(Buffer<DataType, 1> buffer) {
	auto acc = buffer.template get_access<cl::sycl::access::mode::discard_write>();
	for(int i = 0; i < buffer.get_count(); ++i) {
		acc[i] = DataType{};
	}
}

// wrapper for SYCL accessor
template<typename DataType, cl::sycl::access::mode Mode>
class BufferAccessor {
public:
	cl::sycl::accessor<DataType, 1, Mode, cl::sycl::access::target::host_buffer> acc;
	int offset;

	BufferAccessor(Buffer<DataType, 1>& buffer, int offset) : offset(offset), acc(buffer.template get_access<Mode>()) { }
	//~BufferAccessor() {}

	DataType& operator[](int pos) {
		return acc[pos + offset];
	}

	const DataType& operator[](int pos) const {
		return acc[pos + offset];
	}

};

// wrapper for buffer, holds index offset
template<typename DataType>
class BufferWrapper {
private:
	Buffer<DataType, 1> buffer;
	int offset;

public:

	BufferWrapper(Buffer<DataType, 1> buffer, int offset) : buffer(buffer), offset(offset) { }

	template<cl::sycl::access::mode Mode>
	BufferAccessor<DataType, Mode> get_access() {
		BufferAccessor<DataType, Mode> temp = BufferAccessor<DataType, Mode>(buffer, offset);
		return temp;
	}
};

//int test() {
	//Buffer<double, 1> buffer = Buffer<double, 1>(Range<1>(10));
	//BufferWrapper bufferW = BufferWrapper(buffer, -3);
	//auto acc = bufferW.get_access<cl::sycl::access::mode::write>();
	//acc[0] = 1.0;
//}