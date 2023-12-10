#ifndef HDF5_TYPE_MAP
#define HDF5_TYPE_MAP

#include "hdf5.h"
#include <string>
#include <assert.h>

template<typename T>
static hid_t get_hdf5_data_type();

template<> hid_t get_hdf5_data_type<bool>()					{ return H5Tcopy(H5T_NATIVE_INT);		};

template<> hid_t get_hdf5_data_type<char>()					{ return H5Tcopy(H5T_NATIVE_CHAR);		};
template<> hid_t get_hdf5_data_type<unsigned char>()		{ return H5Tcopy(H5T_NATIVE_UCHAR);		};
template<> hid_t get_hdf5_data_type<std::string>()			{
	hid_t datatype = H5Tcopy(H5T_C_S1);
	// Set size to be variable
	auto status = H5Tset_size(datatype, H5T_VARIABLE );
	assert(status >= 0);
	return datatype;     		};
template<> hid_t get_hdf5_data_type<short>()				{ return H5Tcopy(H5T_NATIVE_SHORT);		};
template<> hid_t get_hdf5_data_type<unsigned short>()		{ return H5Tcopy(H5T_NATIVE_USHORT);	};
template<> hid_t get_hdf5_data_type<int>()					{ return H5Tcopy(H5T_NATIVE_INT);		};
template<> hid_t get_hdf5_data_type<unsigned int>()			{ return H5Tcopy(H5T_NATIVE_UINT);		};
template<> hid_t get_hdf5_data_type<long>()					{ return H5Tcopy(H5T_NATIVE_LONG);		};
template<> hid_t get_hdf5_data_type<unsigned long>()		{ return H5Tcopy(H5T_NATIVE_ULONG);		};
template<> hid_t get_hdf5_data_type<long long>()			{ return H5Tcopy(H5T_NATIVE_LLONG);		};
template<> hid_t get_hdf5_data_type<unsigned long long>()	{ return H5Tcopy(H5T_NATIVE_ULLONG);	};

template<> hid_t get_hdf5_data_type<float>()				{ return H5Tcopy(H5T_NATIVE_FLOAT);		};
template<> hid_t get_hdf5_data_type<double>()				{ return H5Tcopy(H5T_NATIVE_DOUBLE);	};
template<> hid_t get_hdf5_data_type<long double>()			{ return H5Tcopy(H5T_NATIVE_LDOUBLE);	};


#endif
