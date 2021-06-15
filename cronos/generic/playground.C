#include <iostream>
#include <queue.H>
#include <matrix.H>

int main() {
	Buffer<double, 1> internal = Buffer<double, 1>(Range<1>(10));
	BufferWrapper<double> buffer = BufferWrapper<double>(internal, -3);
	auto acc = buffer.template get_access<cl::sycl::access::mode::write>();
	//acc[0] = 2.0;

	//auto acc2 = buffer.get_access<cl::sycl::access::mode::read>();
	//double b = acc[0];

	return 0;
}

//int main() {
//
//	auto queue = cl::sycl::queue{};
//
//	const int low = -2;
//	const int high = 14;
//	const int size = high-low;
//	NumMatrix<float, 1> matrix;
//
//	matrix.resize(&low, &high);
//
//	for (int i = low; i < high; ++i) {
//		matrix(i) = i;
//	}
//
//	using SingleDimType = Buffer<float, 1>;
//	std::vector<float> temp(size);
//	for (int i = low; i < high; ++i) {
//		temp[i-low] = i;
//	}
//	SingleDimType buffer(temp.data(), Range<1>(size));
//
//	assert(buffer.get_count() == size);
//
//	//auto acc = buffer.get_access<cl::sycl::access::mode::read>();
//	//for (int i = 0; i < buffer.get_count(); ++i) {
//	//	std::cout << "i: " << i << ", " << acc[i] << std::endl;
//	//}
//
//	bool result = isEqual(buffer, matrix);
//
//	std::cout << "Result of comparison: " << result << std::endl;
//	
//	return 0;
//}

//#include <iostream>
//
//using namespace cl::sycl;
//
//int main(int argc, char* argv[]) {
//
//    /* Create 1024 element arrays for input data and output. */
//    float inputDataA[1024];
//    float inputDataB[1024];
//    float outputData[1024];
//
//    /* Initialize input data with values and output data with zeroes. */
//    for (int i = 0; i < 1024; i++) {
//        inputDataA[i] = (float)i;
//        inputDataB[i] = (float)(1024 - i);
//        outputData[i] = 0.0f;
//    }
//
//    /* Wrap all SYCL structures and function calls with a try-catch block to catch
//     * SYCL exceptions. */
//    try {
//
//        /* Create a default_selector to select a device to execute on. */
//        host_selector mySelector;
//
//        /* Create a queue from the default_selector to create an implicit context
//         * and queue to execute with. */
//        queue myQueue(mySelector);
//
//        /* Create a scope to control data synchronisation of buffer objects. */
//        {
//            /* Create 1 dimensionsal buffers for the input and output data of size
//             * 1024. */
//            buffer<float, 1> inputBufferA(inputDataA, range<1>(1024));
//            buffer<float, 1> inputBufferB(inputDataB, range<1>(1024));
//            buffer<float, 1> outputBuffer(outputData, range<1>(1024));
//
//            /* Submit a command_group to execute from the queue. */
//            myQueue.submit([&](handler& cgh) {
//
//                /* Create accessors for accessing the input and output data within the
//                 * kernel. */
//                auto inputPtrA = inputBufferA.get_access<access::mode::read>(cgh);
//                auto inputPtrB = inputBufferB.get_access<access::mode::read>(cgh);
//                auto outputPtr = outputBuffer.get_access<access::mode::write>(cgh);
//
//                /* Enqueue a kernel called 'vector_add', with a global work size of {
//                 * 16, 8, 8 } and a local work size of { 4, 2, 2 }. */
//                cgh.parallel_for<class vector_add>(
//                    nd_range<3>(range<3>(16, 8, 8), range<3>(4, 2, 2)),
//                    [=](nd_item<3> item) {
//
//                    /* Retreive the linear global id for the current work item. */
//                    size_t idx = item.get_global_linear_id();
//
//                    /* Use the linear global id to add the respective element of each
//                     * input accessor together and assign them to the respective
//                     * element of the output accessor. */
//                    outputPtr[idx] = inputPtrA[idx] + inputPtrB[idx];
//                });
//            });
//        }
//
//    } catch (exception e) {
//
//        /* In the case of an exception being throw, print theerror message and
//         * return 1. */
//        std::cout << e.what();
//        return 1;
//    }
//
//    /* Sum up all the values in the output array. */
//    float sum = 0;
//    for (int i = 0; i < 1024; i++) {
//        sum += outputData[i];
//    }
//
//    /* If the sum is the expected result, return 0, else return 1. */
//    if (sum == (1024.0f * 1024)) {
//        std::cout << "Success!" << std::endl;
//        return 0;
//    } else {
//        std::cout << "Fail: Expected result was 1024.0f, actual result is " << sum
//            << std::endl;
//        return 1;
//    }
//}

//#include <CL/sycl.hpp>
//#include <iostream>
//
//class vector_addition;
//
//int main(int, char**) {
//	cl::sycl::float4 a = { 1.0, 2.0, 3.0, 4.0 };
//	cl::sycl::float4 b = { 4.0, 3.0, 2.0, 1.0 };
//	cl::sycl::float4 c = { 0.0, 0.0, 0.0, 0.0 };
//
//	cl::sycl::host_selector device_selector;
//
//	cl::sycl::queue queue(device_selector);
//	std::cout << "Running on "
//		<< queue.get_device().get_info<cl::sycl::info::device::name>()
//		<< "\n";
//	{
//		cl::sycl::buffer<cl::sycl::float4, 1> a_sycl(&a, cl::sycl::range<1>(1));
//		cl::sycl::buffer<cl::sycl::float4, 1> b_sycl(&b, cl::sycl::range<1>(1));
//		cl::sycl::buffer<cl::sycl::float4, 1> c_sycl(&c, cl::sycl::range<1>(1));
//
//		queue.submit([&](cl::sycl::handler& cgh) {
//			auto a_acc = a_sycl.get_access<cl::sycl::access::mode::read>(cgh);
//			auto b_acc = b_sycl.get_access<cl::sycl::access::mode::read>(cgh);
//			auto c_acc = c_sycl.get_access<cl::sycl::access::mode::discard_write>(cgh);
//
//			cgh.single_task<class vector_addition>([=]() {
//				c_acc[0] = a_acc[0] + b_acc[0];
//			});
//		});
//	}
//	std::cout << "  A { " << a.x() << ", " << a.y() << ", " << a.z() << ", " << a.w() << " }\n"
//		<< "+ B { " << b.x() << ", " << b.y() << ", " << b.z() << ", " << b.w() << " }\n"
//		<< "------------------\n"
//		<< "= C { " << c.x() << ", " << c.y() << ", " << c.z() << ", " << c.w() << " }"
//		<< std::endl;
//
//	return 0;
//}