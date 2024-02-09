#pragma once

#include <CL/sycl.hpp>
#include <celerity/celerity.h>

#define CELERITY_PARALLEL CRONOS_ON

//using Queue = cl::sycl::queue;
using Queue = celerity::distr_queue;

//using DeviceSelector = cl::sycl::host_selector;

template<typename DataType, int Dims>
using CelerityBuffer = celerity::buffer<DataType, Dims>;

template<int Dims>
using CelerityRange = celerity::range<Dims>;

// template<int Dims>
// using Range = cl::sycl::range<Dims>;
