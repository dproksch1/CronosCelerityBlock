#include <stdlib.h>
#include <iostream>
#include "Hdf5File_cbase.H"
#include "Hdf5File_cbase_typemap.H"
#include "buildinfo.H"

// #ifndef H5_NO_NAMESPACE
// using namespace H5;
// #endif
using namespace std;

// -----------------------------------------------------------
// Specializations have to come first (by standard-definition)

template <>
bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, bool AttrData, hid_t my_group) {
  return AddGlobalAttr(AttrName, static_cast<int>(AttrData), my_group);
}

template <>
bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const bool *AttrData, hid_t my_group,
                               int entries) {
  int data[1];
  data[0] = static_cast<int>(AttrData[0]);
  return AddGlobalAttr(AttrName, data, my_group, entries);
}

template <>
bool Hdf5iStream::ReadGlobalAttr(hid_t my_group, const std::string &AttrName, bool &AttrData) {
  int data;
  ReadGlobalAttr(my_group, AttrName, data);
  AttrData = static_cast<bool>(data);
  return true;
}

template <>
bool Hdf5Stream::ChangeGlobalAttr(const std::string &AttrName, bool AttrData, hid_t my_group) {
  return ChangeGlobalAttr(AttrName, static_cast<int>(AttrData), my_group);
}

template <>
bool Hdf5Stream::AddAttrToArrSingle(const std::string &ArrayName, hid_t my_group,
                                    const std::string &AttributeName, bool AttrData) {
  return AddAttrToArrSingle(ArrayName, my_group, AttributeName, static_cast<int>(AttrData));
}

// End Specializations
// -----------------------------------------------------------

// -----------------------------------------------------------
// Begin Hdf5Stream:: Stuff

Hdf5Stream::Hdf5Stream(std::string filename, int NumEntries, int rank, bool use_MPI_IO)  {
  this->filename = filename;
  this->rank = rank;
  this->use_MPI_IO = use_MPI_IO;
  this->NumEntries = NumEntries;

#if (HDF_PARALLEL_IO == CRONOS_ON)
  // Set specific property lists:
  plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
  plist_dset_id = H5Pcreate(H5P_DATASET_XFER);
#else
  this->use_MPI_IO = false;
#endif

  if (!this->use_MPI_IO) {
    plist_file_id = H5P_DEFAULT;
    plist_dset_id = H5P_DEFAULT;
  

    hdf5file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file_id);

    group = H5Gcreate2(hdf5file, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create dataspace
    hid_t info_id = H5Screate(H5S_SCALAR);
    //
    //	// Create Attribute
    //	hid_t info = H5Acreate2(group, "Entries", H5T_NATIVE_INT, info_id,
    //	                        H5P_DEFAULT, H5P_DEFAULT);
    //	// hid_t info = H5Acreate2(group, "Entries", H5T_NATIVE_INT, H5S_SCALAR,
    //	//                         H5P_DEFAULT, H5P_DEFAULT);
    //	// Write Attribute
    //	cout << " Entries set to " << NumEntries << endl;
    //	H5Awrite(info, H5T_NATIVE_INT, &NumEntries);
    //	// Close Attribute
    //	H5Aclose(info);

    // DataSpace infospace( 1, &DimsInfo);
    // datatype.setOrder( H5T_ORDER_LE ); // Little endian
    // Attribute* Info = new Attribute(group->createAttribute("Entries", datatype, infospace));
    // Info->write( PredType::NATIVE_INT, &NumEntries);

    // Indicate new version of hdf5 output:
    // Create Attribute
    hid_t info = H5Acreate2(group, "using_cbase", H5T_NATIVE_INT, info_id, H5P_DEFAULT, H5P_DEFAULT);
    int cbase_val = 1;
    // Write Attribute
    H5Awrite(info, H5T_NATIVE_INT, &cbase_val);
    // Close Attribute
    H5Aclose(info);

    // Close dataspace
    H5Sclose(info_id);

    this->NumEntries = NumEntries;
    this->num = 0;
    this->open = true;

    // Write version and build information
    hid_t version_group = AddGroup("/Data/version");
    AddGlobalAttr("CRONOS_GIT_VERSION", CronosGitVersion(), version_group);
    AddGlobalAttr("CRONOS_GIT_COMMIT", CronosGitCommit(), version_group);
    AddGlobalAttr("APP_GIT_REPO", ApplicationGitRepo(), version_group);
    AddGlobalAttr("APP_GIT_VERSION", ApplicationGitVersion(), version_group);
    AddGlobalAttr("APP_GIT_COMMIT", ApplicationGitCommit(), version_group);
    AddGlobalAttr("BUILD_DATE", BuildDate(), version_group);
printf("close stream_init\n");
    CloseGroup(version_group);
  }
}

void Hdf5Stream::Initialize(MPI_Comm comm) {
  if (use_MPI_IO) {
    // Create file access property list for parallel I/O

    MPI_Info info_mpi = MPI_INFO_NULL;
    // Bug fix for gcc > 6
    //		H5Pset_fapl_mpio(plist_file_id, comm, info_mpi);
    H5Pset_fapl_mpio(plist_file_id, comm, info_mpi);

    // property list for dataset access
    H5Pset_dxpl_mpio(plist_dset_id, H5FD_MPIO_COLLECTIVE);
  } else  {
    plist_file_id = H5P_DEFAULT;
    plist_dset_id = H5P_DEFAULT;
  }

  hdf5file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file_id);

  group = H5Gcreate2(hdf5file, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Create dataspace
  hid_t info_id = H5Screate(H5S_SCALAR);
  //
  //	// Create Attribute
  //	hid_t info = H5Acreate2(group, "Entries", H5T_NATIVE_INT, info_id,
  //	                        H5P_DEFAULT, H5P_DEFAULT);
  //	// hid_t info = H5Acreate2(group, "Entries", H5T_NATIVE_INT, H5S_SCALAR,
  //	//                         H5P_DEFAULT, H5P_DEFAULT);
  //	// Write Attribute
  //	cout << " Entries set to " << NumEntries << endl;
  //	H5Awrite(info, H5T_NATIVE_INT, &NumEntries);
  //	// Close Attribute
  //	H5Aclose(info);

  // DataSpace infospace( 1, &DimsInfo);
  // datatype.setOrder( H5T_ORDER_LE ); // Little endian
  // Attribute* Info = new Attribute(group->createAttribute("Entries", datatype, infospace));
  // Info->write( PredType::NATIVE_INT, &NumEntries);

  // Indicate new version of hdf5 output:
  // Create Attribute
  hid_t info = H5Acreate2(group, "using_cbase", H5T_NATIVE_INT, info_id, H5P_DEFAULT, H5P_DEFAULT);
  int cbase_val = 1;
  // Write Attribute
  H5Awrite(info, H5T_NATIVE_INT, &cbase_val);
  // Close Attribute
  H5Aclose(info);

  // Close dataspace
  H5Sclose(info_id);

  this->num = 0;
  this->open = true;

  // Write version and build information
  hid_t version_group = AddGroup("/Data/version");
  AddGlobalAttr("CRONOS_GIT_VERSION", CronosGitVersion(), version_group);
  AddGlobalAttr("CRONOS_GIT_COMMIT", CronosGitCommit(), version_group);
  AddGlobalAttr("APP_GIT_REPO", ApplicationGitRepo(), version_group);
  AddGlobalAttr("APP_GIT_VERSION", ApplicationGitVersion(), version_group);
  AddGlobalAttr("APP_GIT_COMMIT", ApplicationGitCommit(), version_group);
  AddGlobalAttr("BUILD_DATE", BuildDate(), version_group);

  CloseGroup(version_group);
}

void Hdf5Stream::increase_num() {
  num++;
}

void Hdf5Stream::Reopen() {
  if (!open) {
    hdf5file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, plist_file_id);

    group = H5Gopen2(hdf5file, "Data", H5P_DEFAULT);
  }
  this->open = true;
}

hid_t Hdf5Stream::AddGroup(std::string groupName) {
  //! Add a group to the file
  hid_t new_group = H5Gcreate2(hdf5file, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // Now prepare storage for number of entries
  // Create dataspace
  //	hid_t info_id = H5Screate(H5S_SCALAR);
  //	// Create Attribute
  //	hid_t info = H5Acreate2(new_group, "Entries", H5T_NATIVE_INT, info_id,
  //			H5P_DEFAULT, H5P_DEFAULT);
  //	int num_group(0);
  //	H5Awrite(info, H5T_NATIVE_INT, &num_group);
  //	H5Aclose(info);
  //	H5Sclose(info_id);
  return new_group;
}

bool Hdf5Stream::AddToEntries() { return AddToEntries(group); }
bool Hdf5Stream::AddToEntries(hid_t my_group) {
  hid_t attr_entries;
  int num_entries;
  if (!doesAttrExist(my_group, "Entries")) {
    hid_t info_id = H5Screate(H5S_SCALAR);
    attr_entries =
        H5Acreate2(my_group, "Entries", H5T_NATIVE_INT, info_id, H5P_DEFAULT, H5P_DEFAULT);
    hid_t return_val = H5Sclose(info_id);
    num_entries = 0;
  } else {
    attr_entries = H5Aopen(my_group, "Entries", H5P_DEFAULT);
    H5Aread(attr_entries, H5T_NATIVE_INT, &num_entries);
  }

  // Write attribute
  num_entries++;
  H5Awrite(attr_entries, H5T_NATIVE_INT, &num_entries);
  // Close attribute
  H5Aclose(attr_entries);
  return true;
}

void Hdf5Stream::CloseGroup(hid_t my_group) { H5Gclose(my_group); }

hid_t Hdf5Stream::get_defaultGroup() { return group; }

std::string Hdf5iStream::get_defaultGroupName() { return groupname; }

hid_t Hdf5iStream::get_defaultGroup() { return group; }

hid_t Hdf5iStream::OpenGroup(std::string GroupName) {
  // the default group is opened explicitly in the constructor
  if (GroupName == groupname) {
    return group;
  }

  // Check existence
  if (!doesGroupExist(groupname)) {
    cerr << " There is no group named: " << GroupName << endl;
    exit(-2);
  }
  hid_t my_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  return my_group;
}

void Hdf5iStream::CloseGroup(hid_t my_group) {
  // the default group is closed explicitly in the constructor
  if (my_group != group) {
    H5Gclose(my_group);
  }
}

Hdf5iStream::Hdf5iStream(std::string filename, int rank, bool use_MPI_IO) {
  if (rank == 0) {
    cout << " Opening file: " << filename;
  }

  if (!file_exists(filename)) {
    throw std::invalid_argument("HDF5 File not found: " + filename);
  }

  return_val = 0;

  // Getting the version
  // MinorNum=0;
  // MajorNum=0;
  // ReleaseNum=0;
  // H5::H5Library::getLibVersion(MajorNum, MinorNum, ReleaseNum);

  // int len = filename.size();
  // char* fname = new char[len+1];
  // filename.copy(fname,len);
  // fname[len] = '\0';

  // Open File

#if (HDF_PARALLEL_IO == CRONOS_ON)
  this->use_MPI_IO = use_MPI_IO;

  if (this->use_MPI_IO) {
    // Create file access property list for parallel I/O
    plist_file_id = H5Pcreate(H5P_FILE_ACCESS);
    MPI_Info info_mpi = MPI_INFO_NULL;
    // Bug fix for gcc > 6
    //    H5Pset_fapl_mpio(plist_file_id, comm, info_mpi);
    H5Pset_fapl_mpio(plist_file_id, MPI_COMM_WORLD, info_mpi);

    // property list for dataset access
    plist_dset_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_dset_id, H5FD_MPIO_COLLECTIVE);
  }

#else
  this->use_MPI_IO = false;
#endif

  if (!this->use_MPI_IO) {
    plist_file_id = H5P_DEFAULT;
    plist_dset_id = H5P_DEFAULT;
  }

  // hdf5file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hdf5file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, plist_file_id);

  // Open Group in File
  groupname = "Data";
  group = H5Gopen2(hdf5file, groupname.c_str(), H5P_DEFAULT);

  if (rank == 0) {
    cout << " ...done " << endl;
  }

  //	int numh5;
  //	// Open Attribute "Entries"
  //	hid_t Info = H5Aopen(group, "Entries", H5P_DEFAULT);
  //	// Read Attribute (number of entries)
  //	H5Aread(Info, H5T_NATIVE_INT, &numh5);
  //	H5Aclose(Info);
  //
  //	this->NumEntries = numh5;
}

bool Hdf5Stream::doesAttrExist(hid_t h5group, std::string name) const {
  return (H5Aexists(h5group, name.c_str()) > 0 ? true : false);
}

bool Hdf5Stream::AddDatasetName(const std::string &DatasetName) {
  return AddDatasetName(DatasetName, group);
}

bool Hdf5Stream::AddDatasetName(const std::string &DatasetName, hid_t my_group) {
  int numEntries;
  // Advance number of entries and write corresponding dataset
  hid_t attr_entries;

  if (!doesAttrExist(my_group, "Entries")) {
    hid_t info_id = H5Screate(H5S_SCALAR);
    attr_entries =
        H5Acreate2(my_group, "Entries", H5T_NATIVE_INT, info_id, H5P_DEFAULT, H5P_DEFAULT);
    hid_t return_val = H5Sclose(info_id);
    numEntries = 0;
  } else {
    attr_entries = H5Aopen(my_group, "Entries", H5P_DEFAULT);
    H5Aread(attr_entries, H5T_NATIVE_INT, &numEntries);
  }

  // Now add name of the group
  char cnum[255];
  sprintf(cnum, "%2.2d", numEntries);

  std::string AttrName = "Name_om";
  AttrName += cnum;

  // Set datatype to std::string
  hid_t datatype = H5Tcopy(H5T_C_S1);

  // get length of std::string
  hsize_t StrLen = DatasetName.size();

  // Avoid error for empty std::strings
  if (StrLen == 0) {
    ++StrLen;
  }

  // Set size to length of std::string
  H5Tset_size(datatype, StrLen);

  // Set order to little endian
  H5Tset_order(datatype, H5T_ORDER_LE);

  // Create dataspace for attribute
  hid_t AttrSpace = H5Screate(H5S_SCALAR);

  // Create Attribute
  hid_t info =
      H5Acreate2(my_group, AttrName.c_str(), datatype, AttrSpace, H5P_DEFAULT, H5P_DEFAULT);

  // Write attribute
  hid_t return_val = H5Awrite(info, datatype, DatasetName.c_str());

  // Now increase number of entries
  numEntries++;
  return_val = H5Awrite(attr_entries, H5T_NATIVE_INT, &numEntries);
  return_val = H5Aclose(attr_entries);

  // Close Dataspace
  H5Sclose(AttrSpace);
  // Close Attribute
  H5Aclose(info);
  return true;
}

std::string Hdf5iStream::GetDatasetName(std::string GroupName, int num) {
  hid_t loc_group;

  if (GroupName != "Data") {
    loc_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  } else {
    loc_group = group;
  }

  char cnum[255];
  sprintf(cnum, "%2.2d", num);
  std::string AttrName = "Name_om";
  AttrName += cnum;

  if (!doesAttrExist(AttrName, loc_group)) {
    if (GroupName != "Data") {
      // Close group
      H5Gclose(loc_group);
    }
    return "";
  }

  // Open Attribute
  hid_t Info = H5Aopen(loc_group, AttrName.c_str(), H5P_DEFAULT);
  hid_t ftype = H5Aget_type(Info);
  hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);

  htri_t size_var;
  std::string DatasetName;

  if ((size_var = H5Tis_variable_str(ftype)) > 0) {
    char string_attr[256];
    H5Aread(Info, type, &string_attr);
    DatasetName = string_attr;
  } else {
    char string_out[255];

    size_t sdim = H5Tget_size(ftype);
    sdim++; /* Make room for null terminator */
    hid_t space = H5Aget_space(Info);
    hsize_t dims[1];
    hid_t ndims = H5Sget_simple_extent_dims(space, dims, NULL);

    hid_t memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, sdim);

    H5Aread(Info, memtype, string_out);
    DatasetName = string_out;
  }
  H5Aclose(Info);

  if (GroupName != "Data") {
    // Close group
    H5Gclose(loc_group);
  }

  return DatasetName;
}

std::string Hdf5iStream::GetDatasetName(int num) { return GetDatasetName("Data", num); }

// std::string Hdf5iStream::GetDatasetName(int num) {
//
//	char cnum[255];
//	sprintf(cnum,"%2.2d",num);
//	const std::string &AttrName = "Name_om";
//	AttrName += cnum;
//
//
//	// Open Attribute
//	hid_t Info = H5Aopen(group, AttrName.c_str(), H5P_DEFAULT);
//	hid_t ftype = H5Aget_type(Info);
//	hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
//
//	htri_t size_var;
//	std::string DatasetName;
//
//	if((size_var = H5Tis_variable_str(ftype)) == 1) {
//		char *string_attr;
//		H5Aread(Info, type, &string_attr);
//		DatasetName = string_attr;
//	} else {
//		char string_out[255];
//		H5Aread(Info, type, string_out);
//		DatasetName = string_out;
//	}
//	H5Aclose(Info);
//
//	return DatasetName;
//}

#if (HDF_PARALLEL_IO == CRONOS_ON)
template <typename T, typename>
bool Hdf5Stream::WriteNumArray_withMPI_IO(const std::string &ArrayName, const NumArray<T> &data,
                                          bool with_data) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

  */
  int mx = data.getLength();

  num += 1;

  int DIM = 1;

  // Choose float, little endian of size 1
  hid_t datatype = get_hdf5_data_type<T>();
  return_val = H5Tset_order(datatype, H5T_ORDER_LE);

  // Create dataspace
  hsize_t DimsData = mx;
  hid_t dataspace = H5Screate_simple(DIM, &DimsData, NULL);

  // Create dataset
  hid_t dataset = H5Dcreate(group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
  H5Sclose(dataspace);

  hid_t memspace = H5Screate_simple(DIM, &DimsData, NULL);
  if (!with_data) {
    H5Sselect_none(memspace);
  }

  hid_t groupspace = H5Dget_space(dataset);

  hsize_t offset[1] = {0};
  H5Sselect_hyperslab(groupspace, H5S_SELECT_SET, offset, NULL, &DimsData, NULL);
  if (!with_data) {
    H5Sselect_none(groupspace);
  }

  // Prepare mpi-stuff
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // Write data
  //	if(with_data) {
  return_val = H5Dwrite(dataset, datatype, memspace, groupspace, plist_id, data);
  //	}

  H5Dclose(dataset);
  H5Pclose(plist_id);
  H5Sclose(memspace);
  H5Sclose(groupspace);

  return true;
}

template <typename T, typename>
bool Hdf5Stream::WriteNumArray_withMPI_IO2(const std::string &ArrayName, const NumArray<T> &data,
                                           bool with_data) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

  */
  int mx = data.getLength();

  num += 1;

  int DIM = 1;

  // Choose float, little endian of size 1
  hid_t datatype = get_hdf5_data_type<T>();
  return_val = H5Tset_order(datatype, H5T_ORDER_LE);

  // Create dataspace
  hsize_t DimsData = mx;
  hid_t dataspace = H5Screate_simple(DIM, &DimsData, NULL);

  hid_t memspace = H5Screate_simple(DIM, &DimsData, NULL);
  if (!with_data) {
    H5Sselect_none(memspace);
  }

  hsize_t offset[1] = {0};
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, &DimsData, NULL);
  if (!with_data) {
    H5Sselect_none(dataspace);
  }

  // Create dataset
  hid_t dataset = H5Dcreate2(group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);

  // Prepare mpi-stuff
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  // Write data
  //	if(with_data) {
  return_val = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, plist_id, data);
  //	}
  H5Pclose(plist_id);

  cout << " Writing float " << endl;

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(memspace);
  return true;
}

#endif

int Hdf5iStream::GetDatasetDimension(std::string DataSetName) {
  hid_t dataset = H5Dopen2(group, DataSetName.c_str(), H5P_DEFAULT);

  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset);
  int dimhdf = H5Sget_simple_extent_ndims(dataspace);

  H5Dclose(dataset);

  return dimhdf;
}

NumArray<int> Hdf5iStream::GetDatasetExtent(std::string DataSetName) {
  //! Get extent of each dimension
  hid_t dataset = H5Dopen2(group, DataSetName.c_str(), H5P_DEFAULT);

  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset);

  int ndims = H5Sget_simple_extent_ndims(dataspace);

  std::vector<hsize_t> dims_out(ndims);

  ndims = H5Sget_simple_extent_dims(dataspace, dims_out.data(), NULL);

  NumArray<int> Nx(ndims);

  cout << endl << " my dims: " << ndims << " " << endl << endl;

  for (int idim = 0; idim < ndims; ++idim) {
    Nx[idim] = dims_out[ndims - (idim + 1)];
    // Nx[idim] = dims_out[idim];
  }

  H5Dclose(dataset);

  return Nx;
}

// int Hdf5iStream::ReadDatasetAttr(std::string DataSetName, const std::string &AttrName)
// {
// 	// Open the dataset
// 	hid_t dataset = H5Dopen2(group, DataSetName.c_str(), H5P_DEFAULT);

// 	hid_t attr = H5Aopen(dataset, AttrName.c_str(), H5P_DEFAULT);
// 	int value;
// 	H5Aread(attr, H5T_NATIVE_INT, &value);
// 	return value;

// }

bool Hdf5iStream::ReadDatasetGrid(std::string DataSetName, NumArray<float> &xPos,
                                  NumArray<float> &yPos, NumArray<float> &zPos) {
  // Open the dataset
  hid_t dataset = H5Dopen2(group, DataSetName.c_str(), H5P_DEFAULT);

  float xb[3];
  hid_t Origin = H5Aopen(dataset, "origin", H5P_DEFAULT);
  H5Aread(Origin, H5T_NATIVE_FLOAT, xb);
  H5Aclose(Origin);

  float dx[3];
  hid_t Delta = H5Aopen(dataset, "delta", H5P_DEFAULT);
  H5Aread(Delta, H5T_NATIVE_FLOAT, dx);
  H5Aclose(Delta);

  hid_t dataspace = H5Dget_space(dataset);
  int ndims = H5Sget_simple_extent_ndims(dataspace);
  std::vector<hsize_t> dims_out(ndims);
  ndims = H5Sget_simple_extent_dims(dataspace, dims_out.data(), NULL);

  // Now build the data-arrays:
  xPos.resize(dims_out[2]);
  yPos.resize(dims_out[1]);
  zPos.resize(dims_out[0]);

  for (unsigned int ix = 0; ix < dims_out[2]; ++ix) {
    xPos[ix] = xb[0] + dx[0] * ix;
  }
  for (unsigned int iy = 0; iy < dims_out[1]; ++iy) {
    yPos[iy] = xb[1] + dx[1] * iy;
  }
  for (unsigned int iz = 0; iz < dims_out[0]; ++iz) {
    zPos[iz] = xb[2] + dx[2] * iz;
  }

  // Close the dataset:
  H5Dclose(dataset);

  return true;
}

void Hdf5iStream::getSize(std::string GroupName, const std::string &ArrayName, int Nx[], int DIM) {
  // Open group
  hid_t loc_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);

  hid_t dataset = H5Dopen2(loc_group, ArrayName.c_str(), H5P_DEFAULT);
  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset);

  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }

  std::vector<hsize_t> dims_out(DIM);
  H5Sget_simple_extent_dims(dataspace, dims_out.data(), NULL);
  for (int i = 0; i < DIM; ++i) {
    Nx[DIM - i - 1] = dims_out[i];
  }
  H5Dclose(dataset);

  // Close group
  H5Gclose(loc_group);
}

void Hdf5iStream::getSize(const std::string &ArrayName, int mx[], int DIM) {
  hid_t dataset = H5Dopen2(group, ArrayName.c_str(), H5P_DEFAULT);
  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset);

  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }

  std::vector<hsize_t> dims_out(DIM);
  H5Sget_simple_extent_dims(dataspace, dims_out.data(), NULL);
  for (int i = 0; i < DIM; ++i) {
    mx[DIM - i - 1] = dims_out[i];
  }
  H5Dclose(dataset);
}

#if (HDF_PARALLEL_IO == CRONOS_ON)

template <typename T, typename>
bool Hdf5Stream::Write3DMatrix_withMPI_IO(const std::string &ArrayName, const NumMatrix<T, 3> &data,
                                          const NumArray<int> &mx_global,
                                          const NumArray<int> &mx_local,
                                          const NumArray<int> &rank_shift,
                                          // NumArray<int> &rank_pos,
                                          const NumArray<float> &xb, const NumArray<float> &dx) {
  return Write3DMatrix_withMPI_IO(ArrayName, data, mx_global, mx_local, rank_shift, xb, dx, group, true);
}

template <typename T, typename>
bool Hdf5Stream::Write3DMatrix_withMPI_IO(const std::string &ArrayName, const NumMatrix<T, 3> &data,
                                          const NumArray<int> &mx_global,
                                          const NumArray<int> &mx_local,
                                          const NumArray<int> &rank_shift,
                                          // NumArray<int> &rank_pos,
                                          const NumArray<float> &xb, const NumArray<float> &dx,
                                          hid_t my_group, int q_index, bool with_opendxinfo) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file. MPI
     parallel form of the routine that writes all data to a single file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

     rim can be computed from mx_local.

  */

  // Determine rim of data
  int rim = 0;
  if (data.getLow(0) < 0) {
    rim = -data.getLow(0);
  }

  int mx[3];
  mx[0] = data.getHigh(2) - data.getLow(2) + 1;
  mx[1] = data.getHigh(1) - data.getLow(2) + 1;
  mx[2] = data.getHigh(0) - data.getLow(2) + 1;

  AddDatasetName(ArrayName, my_group);
  num += 1;

  int DIM = 3;
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  hsize_t DimsData[DIM];
  DimsData[0] = mx_global[2] + 1 + 2 * rim;
  DimsData[1] = mx_global[1] + 1 + 2 * rim;
  DimsData[2] = mx_global[0] + 1 + 2 * rim;

  hid_t dataspace = H5Screate_simple(DIM, DimsData, NULL);

  // Supplying additional attributes for opendx input

  hid_t datatypefloat = get_hdf5_data_type<float>();
  H5Tset_order(datatypefloat, H5T_ORDER_LE);

  // Create dataspace for attribute
  hsize_t DimsAttr = 3;
  hid_t attrspace = H5Screate_simple(1, &DimsAttr, NULL);
  float Origin[3];
  float Delta[3];
  if (with_opendxinfo) {
    for (int q = 0; q < 3; ++q) {
      Origin[q] = xb[q];
      Delta[q] = dx[q];
    }
  }

  // Create dataset
  hid_t dataset_id = H5Dcreate2(my_group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
  if (with_opendxinfo) {
    // Create attributes
    hid_t origin =
        H5Acreate2(dataset_id, "origin", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    hid_t delta =
        H5Acreate2(dataset_id, "delta", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);

    // Write attributes
    H5Awrite(origin, datatypefloat, &Origin);
    H5Awrite(delta, datatypefloat, &Delta);

    // Close attributes
    H5Aclose(origin);
    H5Aclose(delta);
    H5Sclose(attrspace);
  }

  if (q_index > 0) {
    // Create dataspace
    hid_t info_id = H5Screate(H5S_SCALAR);
    // Create Attribute
    hid_t info =
        H5Acreate2(dataset_id, "q_index", H5T_NATIVE_INT, info_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(info, H5T_NATIVE_INT, &q_index);
    H5Aclose(info);
    H5Sclose(info_id);
  }

  // Now make distinguis local from global data via hyperslabs:
  hsize_t sizeLocal[DIM];
  sizeLocal[0] = mx_local[2] + 1 + 2 * rim;
  sizeLocal[1] = mx_local[1] + 1 + 2 * rim;
  sizeLocal[2] = mx_local[0] + 1 + 2 * rim;

  hsize_t offset[DIM];
  // offset[0] = rank_pos[2]*rank_shift[2];
  // offset[1] = rank_pos[1]*rank_shift[1];
  // offset[2] = rank_pos[0]*rank_shift[0];
  offset[0] = rank_shift[2];
  offset[1] = rank_shift[1];
  offset[2] = rank_shift[0];

  hid_t dataspaceLocal = H5Screate_simple(DIM, sizeLocal, NULL);

  // Select local data as hyperslab in dataset
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, sizeLocal, NULL);

  // Write data
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dataset_id, datatype, dataspaceLocal, dataspace, plist_id, data);
  H5Pclose(plist_id);

  H5Sclose(dataspaceLocal);
  H5Sclose(dataspace);

  // Close the dataset:
  H5Dclose(dataset_id);

  //	if(rank==0) {
  //		AddToEntries(my_group);
  //	}

  return true;
}

template <typename T, typename>
bool Hdf5Stream::Write2DMatrix_withMPI_IO(const std::string &ArrayName, const NumMatrix<T, 2> &data,
                                          const NumArray<int> &mx_global,
                                          const NumArray<int> &mx_local,
                                          const NumArray<int> &rank_shift,
                                          const NumArray<float> &xb, const NumArray<float> &dx) {
  return Write2DMatrix_withMPI_IO(ArrayName, data, mx_global, mx_local, rank_shift, xb, dx, group,
                                  true);
}

template <typename T, typename>
bool Hdf5Stream::Write2DMatrix_withMPI_IO(const std::string &ArrayName, const NumMatrix<T, 2> &data,
                                          const NumArray<int> &mx_global,
                                          const NumArray<int> &mx_local,
                                          const NumArray<int> &rank_shift,
                                          const NumArray<float> &xb, const NumArray<float> &dx,
                                          hid_t my_group, bool with_opendxinfo) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file. MPI
     parallel form of the routine that writes all data to a single file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

     rim can be computed from mx_local.

  */

  // Determine rim of data
  int rim = 0;
  if (data.getLow(0) < 0) {
    rim = -data.getLow(0);
  }

  AddDatasetName(ArrayName, my_group);
  num += 1;

  int DIM = 2;
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  hsize_t DimsData[DIM];
  DimsData[0] = mx_global[1] + 1 + 2 * rim;
  DimsData[1] = mx_global[0] + 1 + 2 * rim;

  hid_t dataspace = H5Screate_simple(DIM, DimsData, NULL);

  // Supplying additional attributes for opendx input

  hid_t datatypefloat = get_hdf5_data_type<float>();
  H5Tset_order(datatypefloat, H5T_ORDER_LE);

  // Create dataspace for attribute
  hsize_t DimsAttr = 2;
  hid_t attrspace = H5Screate_simple(1, &DimsAttr, NULL);
  float Origin[2];
  float Delta[2];
  if (with_opendxinfo) {
    for (int q = 0; q < 2; ++q) {
      Origin[q] = xb[q];
      Delta[q] = dx[q];
    }
  }

  // Create dataset
  hid_t dataset_id = H5Dcreate2(my_group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
  if (with_opendxinfo) {
    // Create attributes
    hid_t origin =
        H5Acreate2(dataset_id, "origin", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    hid_t delta =
        H5Acreate2(dataset_id, "delta", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);

    // Write attributes
    H5Awrite(origin, datatypefloat, &Origin);
    H5Awrite(delta, datatypefloat, &Delta);

    // Close attributes
    H5Aclose(origin);
    H5Aclose(delta);
    H5Sclose(attrspace);
  }

  // Now make distinguis local from global data via hyperslabs:
  hsize_t sizeLocal[DIM];
  sizeLocal[0] = mx_local[1] + 1 + 2 * rim;
  sizeLocal[1] = mx_local[0] + 1 + 2 * rim;

  hsize_t offset[DIM];
  offset[0] = rank_shift[1];
  offset[1] = rank_shift[0];

  hid_t dataspaceLocal = H5Screate_simple(DIM, sizeLocal, NULL);

  // Select local data as hyperslab in dataset
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, sizeLocal, NULL);

  // Write data
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dataset_id, datatype, dataspaceLocal, dataspace, plist_id, data);
  H5Pclose(plist_id);

  H5Sclose(dataspaceLocal);
  H5Sclose(dataspace);

  // Close the dataset:
  H5Dclose(dataset_id);

  //	if(rank==0) {
  //		AddToEntries(my_group);
  //	}

  return true;
}

#endif

// bool Hdf5iStream::Read3DMatrix(const std::string &ArrayName, NumMatrix<float,3> &data)
//{
//	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.
//
//	   Remarks:
//
//	   On reading the hdf data the swapped dimensions have to be taken
//	   into account
//
//	*/
//
//	int DIM = 3;
//	int lbound[3], ubound[3];
//
//	hid_t dataset = H5Dopen2(group, ArrayName.c_str(), H5P_DEFAULT);
//
//	// Get dataspace handler
//	hid_t dataspace = H5Dget_space(dataset);
//
//	int  dimhdf = H5Sget_simple_extent_ndims(dataspace);
//	if(dimhdf != DIM){
//		cerr << " Wrong dimensionality of input data: " << endl;
//		cerr << dimhdf << " " << DIM << endl;
//		exit(-2);
//	}
//	hsize_t dims_out[DIM];
//	int ndims = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
//
//	if(ndims != DIM) {
//		cerr << " Wrong number of dimensions " << ndims << " - " << DIM << endl;
//		exit(-22);
//	}
//
//	for(int i=0;i<3;i++){
//		lbound[i]=0;
//		ubound[i]=int(dims_out[DIM-(i+1)])-1;
//	}
//	data.resize(lbound,ubound);
//
//	H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
//	// Close the dataset:
//	H5Dclose(dataset);
//
//	return true;
//}

bool Hdf5Stream::close() {
  //	// Reopen attribute
  //	hid_t Info = H5Aopen(group, "Entries", H5P_DEFAULT);
  //	// Write attribute
  //	H5Awrite(Info, H5T_NATIVE_INT, &num);
  //	// Close attribute
  //	H5Aclose(Info);

  // Close group and file
  H5Gclose(group);
  H5Fclose(hdf5file);
  this->open = false;
  return true;
}

bool Hdf5iStream::close() {
  //  delete group;
  // H5Gclose(group);
  // H5Fclose(hdf5file);
  return true;
}

Hdf5Stream::~Hdf5Stream() {
  if (this->use_MPI_IO) {
    H5Pclose( plist_file_id );
    H5Pclose( plist_dset_id );
  }
  if (open) {
    close();
  }
}

Hdf5iStream::~Hdf5iStream() {
  //  delete group and file;
  H5Gclose(group);
  H5Fclose(hdf5file);
}

bool Hdf5iStream::doesAttrExist(std::string name, hid_t group) const {
  return (H5Aexists(group, name.c_str()) > 0 ? true : false);
}

bool Hdf5iStream::doesAttrExist(std::string name) const { return doesAttrExist(name, group); }

bool Hdf5iStream::doesGroupExist(std::string name) const {
  return (H5Lexists(hdf5file, name.c_str(), H5P_DEFAULT) > 0 ? true : false);
}

bool Hdf5iStream::doesSubGroupExist(std::string name) const {
  return (H5Lexists(group, name.c_str(), H5P_DEFAULT) > 0 ? true : false);
}

bool Hdf5iStream::doesDsetExist(std::string DsetName) const {
  htri_t dataset_status = H5Lexists(group, DsetName.c_str(), H5P_DEFAULT);
  return (dataset_status > 0);
}

bool Hdf5iStream::doesDsetExist(std::string GroupName, std::string DsetName) const {
  hid_t loc_group;
  if (GroupName != "Data") {
    loc_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  } else {
    loc_group = group;
  }
  htri_t dataset_status = H5Lexists(loc_group, DsetName.c_str(), H5P_DEFAULT);
  return (dataset_status > 0);
}

// ---------------------------------------------------------
// Beginn Hdf5Stream Templates

template <typename T, typename>
bool Hdf5Stream::ChangeGlobalAttr(const std::string &AttrName, T AttrData) {
  return ChangeGlobalAttr(AttrName, AttrData, group);
}

template <typename T, typename>
bool Hdf5Stream::ChangeGlobalAttr(const std::string &AttrName, T AttrData, hid_t my_group) {
  // Open Attribute
  hid_t info = H5Aopen(my_group, AttrName.c_str(), H5P_DEFAULT);
  // Rewrite Attribute
  H5Awrite(info, get_hdf5_data_type<T>(), &AttrData);
  H5Aclose(info);
  return true;
}

template <typename T, typename>
bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, T AttrData) {
  return AddGlobalAttr(AttrName, AttrData, group);
}

template <typename T, typename>
bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, T AttrData, hid_t my_group) {
  // little endian of size 1
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  // Create dataspace
  hid_t AttrSpace = H5Screate(H5S_SCALAR);

  // Create Attribute
  hid_t info =
      H5Acreate2(my_group, AttrName.c_str(), datatype, AttrSpace, H5P_DEFAULT, H5P_DEFAULT);
  // Write Attribute
  if constexpr (std::is_same<T, std::string>::value) {
      const char* temp = AttrData.c_str();
      H5Awrite(info, datatype, &temp);
  } else {
      H5Awrite(info, datatype, &AttrData);
  }

  // Close Attribute
  H5Aclose(info);
  // Close Dataspace
  H5Sclose(AttrSpace);
  return true;
}

template <typename T, typename>
bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const T *AttrData, int num) {
  return AddGlobalAttr(AttrName, AttrData, group, num);
}

template <typename T, typename>
bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const T *AttrData, hid_t my_group,
                               int entries) {
  // little endian of size num
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  // Create dataspace
  hid_t AttrSpace = H5Screate(H5S_SIMPLE);
  hsize_t dimAttr = entries;
  return_val = H5Sset_extent_simple(AttrSpace, 1, &dimAttr, NULL);

  // Create Attribute
  hid_t info =
      H5Acreate2(my_group, AttrName.c_str(), datatype, AttrSpace, H5P_DEFAULT, H5P_DEFAULT);
  // Write Attribute
  H5Awrite(info, datatype, AttrData);

  // Close Attribute
  H5Aclose(info);
  // Close Dataspace
  H5Sclose(AttrSpace);
  return true;
}

bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, char const *AttrData) {
  return AddGlobalAttr(AttrName, std::string(AttrData));
}

template <typename T, typename>
bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName, const std::string &AttributeName,
                                     T AttributeData) {
  return AddAttributeToArray(ArrayName, AttributeName, AttributeData, group);
}

template <typename T, typename>
bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName, const std::string &AttributeName,
                                     T AttributeData, hid_t my_group) {
  return AddAttrToArrSingle(ArrayName, my_group, AttributeName, AttributeData);
}

template <typename T, typename>
bool Hdf5Stream::AddAttrToArrSingle(const std::string &ArrayName, hid_t my_group,
                                    const std::string &AttributeName, T AttributeData) {
  /*******************************************************
   *  Routine to add Attribute to Any written Array Data *
   *  the corresponding array is identified by its name  *
   ******************************************************/

  // Reopen written dataset:
  hid_t dataset = H5Dopen2(my_group, ArrayName.c_str(), H5P_DEFAULT);

  hid_t AttrSpace = H5Screate(H5S_SCALAR);

  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  // Create Attribute
  hid_t info =
      H5Acreate2(dataset, AttributeName.c_str(), datatype, AttrSpace, H5P_DEFAULT, H5P_DEFAULT);
  // Write Attribute
  H5Awrite(info, datatype, &AttributeData);

  // Close everyting
  H5Sclose(AttrSpace);
  H5Aclose(info);
  H5Dclose(dataset);

  return true;
}

template <typename T, typename>
bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName, const std::string &AttributeName,
                                     const T *AttributeData, const hsize_t entries) {
  return AddAttributeToArray(ArrayName, AttributeName, AttributeData, entries, group);
}

template <typename T, typename>
bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName, const std::string &AttributeName,
                                     const T *AttributeData, const hsize_t entries,
                                     hid_t my_group) {
  /*******************************************************
   *  Routine to add Attribute to Any written Array Data *
   *  the corresponding array is identified by its name  *
   ******************************************************/

  // Reopen written dataset:
  hid_t dataset = H5Dopen2(my_group, ArrayName.c_str(), H5P_DEFAULT);

  hid_t AttrSpace = H5Screate_simple(1, &entries, NULL);

  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  // Create Attribute
  hid_t info =
      H5Acreate2(dataset, AttributeName.c_str(), datatype, AttrSpace, H5P_DEFAULT, H5P_DEFAULT);
  // Write Attribute
  H5Awrite(info, datatype, AttributeData);

  // Close everyting
  H5Sclose(AttrSpace);
  H5Aclose(info);
  H5Dclose(dataset);

  return true;
}

template <typename T, typename>
bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName, const NumMatrix<T, 1> &data,
                               double Origin, double Delta, int numin) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

  */
  int mx = data.getHigh(0) - data.getLow(0) + 1;

  num += 1;

  int DIM = 1;

  // Choose double, little endian of size 1
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  // Create dataspace
  hsize_t DimsData = mx;
  hid_t dataspace = H5Screate_simple(DIM, &DimsData, NULL);

  // Supplying additional attributes for opendx input

  // Datatype: double, little endian of size 1
  hid_t datatypefloat = get_hdf5_data_type<double>();
  H5Tset_order(datatypefloat, H5T_ORDER_LE);

  // Create dataspace for attribute
  hid_t attrspace = H5Screate(H5S_SCALAR);

  // Create dataset
  hid_t dataset = H5Dcreate2(group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
  // Create attributes
  hid_t origin = H5Acreate2(dataset, "origin", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  hid_t delta = H5Acreate2(dataset, "delta", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);

  // Write attributes
  H5Awrite(origin, datatype, &Origin);
  H5Awrite(delta, datatype, &Delta);

  H5Aclose(origin);
  H5Aclose(delta);

  // Write data only if rank 0

#if (HDF_PARALLEL_IO == CRONOS_ON)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
#endif

    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  // Close dataset
  H5Dclose(dataset);

  // Close dataspace
  H5Sclose(dataspace);

  return true;
}

template <typename T, typename>
bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName, const NumMatrix<T, 1> &data) {
  return Write1DMatrix(ArrayName, data, group);
}

template <typename T, typename>
bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName, const NumMatrix<T, 1> &data,
                               hid_t my_group) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

  */
  int mx = data.getHigh(0) - data.getLow(0) + 1;

  num += 1;

  int DIM = 1;

  // Choose float, little endian of size 1
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  // Create dataspace
  hsize_t DimsData = mx;
  hid_t dataspace = H5Screate_simple(DIM, &DimsData, NULL);

  // Create dataset
  hid_t dataset = H5Dcreate2(my_group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);

  // Write data only if rank 0

#if (HDF_PARALLEL_IO == CRONOS_ON)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
#endif

    H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  // Close dataset
  H5Dclose(dataset);

  // Close dataspace
  H5Sclose(dataspace);

  return true;
}

template <typename T, typename>
bool Hdf5Stream::WriteNumArray(const std::string &ArrayName, const NumArray<T> &data) {
  return WriteNumArray(ArrayName, data, group);
}

template <typename T, typename>
bool Hdf5Stream::WriteNumArray(const std::string &ArrayName, const NumArray<T> &data,
                               hid_t my_group) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

  */
  int mx = data.getLength();

  num += 1;

  int DIM = 1;

  // Choose float, little endian of size 1
  hid_t datatype = get_hdf5_data_type<T>();
  return_val = H5Tset_order(datatype, H5T_ORDER_LE);

  // Create dataspace
  hsize_t DimsData = mx;
  hid_t dataspace = H5Screate_simple(DIM, &DimsData, NULL);

  // Create dataset
  hid_t dataset = H5Dcreate2(my_group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);

  // Write data
  return_val = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dataset);
  H5Sclose(dataspace);
  return true;
}

template <typename T, typename>
bool Hdf5Stream::Write3DMatrix(const std::string &ArrayName, const NumMatrix<T, 3> &data,
                               const double *xb, const double *dx) {
  return Write3DMatrix(ArrayName, data, xb, dx, group, true);
}

template <typename T, typename>
bool Hdf5Stream::Write3DMatrix(const std::string &ArrayName, const NumMatrix<T, 3> &data) {
  double dummy[3];
  return Write3DMatrix(ArrayName, data, dummy, dummy, group, false);
}

template <typename T, typename>
bool Hdf5Stream::Write3DMatrix(const std::string &ArrayName, const NumMatrix<T, 3> &data,
                               const double *xb, const double *dx, hid_t my_group, int q_index,
                               bool with_opendxinfo) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

  */
  int mx[3];
  mx[0] = data.getHigh(2) - data.getLow(2) + 1;
  mx[1] = data.getHigh(1) - data.getLow(1) + 1;
  mx[2] = data.getHigh(0) - data.getLow(0) + 1;

  // Add name of dataset to group
  if (my_group != group) {
    AddDatasetName(ArrayName, my_group);
  }
  // Add name globally (deprecated)
  AddDatasetName(ArrayName);
  num += 1;

  const int DIM = 3;
  // Choose double, little endian of size 1
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  hsize_t DimsData[DIM];
  for (int q = 0; q < DIM; ++q) {
    DimsData[q] = mx[q];
  }
  hid_t dataspace = H5Screate_simple(DIM, DimsData, NULL);
  // Supplying additional attributes for opendx input
  // Datatype: double, little endian of size 1
  hid_t datatypefloat = get_hdf5_data_type<double>();
  H5Tset_order(datatypefloat, H5T_ORDER_LE);

  // Create dataspace for attribute
  hsize_t DimsAttr = 3;
  hid_t attrspace = H5Screate_simple(1, &DimsAttr, NULL);

  double Origin[3];
  double Delta[3];
  if (with_opendxinfo) {
    for (int q = 0; q < 3; ++q) {
      Origin[q] = float(xb[q]);
      Delta[q] = float(dx[q]);
    }
  }

  // Create dataset
  hid_t dataset = H5Dcreate2(my_group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);

  if (with_opendxinfo) {
    // Create attributes
    hid_t origin =
        H5Acreate2(dataset, "origin", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    hid_t delta = H5Acreate2(dataset, "delta", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);

    // Write attributes
    H5Awrite(origin, datatypefloat, &Origin);
    H5Awrite(delta, datatypefloat, &Delta);
    // Close attributes
    H5Aclose(origin);
    H5Aclose(delta);
    H5Sclose(attrspace);
  }

  if (q_index > 0) {
    // Create dataspace
    hid_t info_id = H5Screate(H5S_SCALAR);
    // Create Attribute
    hid_t info = H5Acreate2(dataset, "q_index", H5T_NATIVE_INT, info_id, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(info, H5T_NATIVE_INT, &q_index);
    H5Aclose(info);
    H5Sclose(info_id);
  }

  // Write data
  H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);

  H5Sclose(dataspace);

  return true;
}

template <typename T, typename>
bool Hdf5Stream::WriteNDArray(const std::string &ArrayName, const T *data, const int mx[],
                              int dim) {
  // Add name globally (deprecated)
  AddDatasetName(ArrayName);

  num += 1;
  //	if(num > NumEntries){
  //		cerr << " Error (Hdf5Stream::WriteNDArray): more datasets written than allowed for
  // by construcor " << endl; 		return false;
  //	}

  // Set datatype to float
  hid_t datatype = get_hdf5_data_type<T>();
  // Use little endian:
  H5Tset_order(datatype, H5T_ORDER_LE);

  std::vector<hsize_t> DimsData(dim);
  for (int dir = 0; dir < dim; ++dir) {
    DimsData[dir] = mx[dir];
  }

  // Create a dataspace of arbitrary dimension
  hid_t dataspace = H5Screate_simple(dim, DimsData.data(), NULL);

  // Create dataset
  hid_t dataset = H5Dcreate2(group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
  // Write data
  H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);

  // Close dataspace
  H5Sclose(dataspace);

  return true;
}

template <typename T, typename>
bool Hdf5Stream::Write3DVecMatrix(const std::string &ArrayName, const NumMatrix<T, 3> &data_x,
                                  const NumMatrix<T, 3> &data_y, const NumMatrix<T, 3> &data_z) {
  /* Routine to write a Matrix of vectorial values to the file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

  */
  int ib[3];
  ib[0] = data_x.getLow(2);
  ib[1] = data_x.getLow(1);
  ib[2] = data_x.getLow(0);
  int mx[3];
  mx[0] = data_x.getHigh(2) - data_x.getLow(2) + 1;
  mx[1] = data_x.getHigh(1) - data_x.getLow(1) + 1;
  mx[2] = data_x.getHigh(0) - data_x.getLow(0) + 1;

  AddDatasetName(ArrayName);

  num += 1;

  const int DIM = 3;

  /*
* Define array datatype for the data in the file.
- 1 component
- 3 entries
*/
  hsize_t ArrayExtent[1];
  ArrayExtent[0] = 3;
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_size(datatype, 3);
  H5Tset_order(datatype, H5T_ORDER_LE);

  cerr << " Has to be tested " << endl;
  exit(2);

  // ArrayType datatype( PredType::NATIVE_FLOAT, 1, ArrayExtent);
  // datatype.setOrder( H5T_ORDER_LE );

  // FloatType datatype( PredType::NATIVE_DOUBLE );
  // datatype.setOrder( H5T_ORDER_LE );

  hsize_t DimsData[DIM];
  for (int q = 0; q < DIM; ++q) {
    DimsData[q] = mx[q];
  }
  hid_t dataspace = H5Screate_simple(DIM, DimsData, NULL);

  //float data[mx[0]][mx[1]][mx[2]][3];
  std::vector<std::vector<std::vector<std::vector<float>>>> data(mx[0]);
  for (int i = 0; i < mx[0]; i++) {
	  data[i].resize(mx[1]);
	  for (int j = 0; j < mx[1]; j++) {
		  data[i][j].resize(mx[2]);
		  for (int k = 0; k < mx[2]; k++) {
			  data[i][j][k].resize(3);
		  }
	  }
  }
  
  for (int i = 0; i < mx[0]; i++) {
    for (int j = 0; j < mx[1]; j++) {
      for (int k = 0; k < mx[2]; k++) {
        data[i][j][k][0] = data_x(k + ib[2], j + ib[1], i + ib[0]);
        data[i][j][k][1] = data_y(k + ib[2], j + ib[1], i + ib[0]);
        data[i][j][k][2] = data_z(k + ib[2], j + ib[1], i + ib[0]);
      }
    }
  }

  // Create dataset
  hid_t dataset = H5Dcreate2(group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);

  std::vector<float> tempData((double)mx[0]*mx[1]*mx[2]*3);
  long cur = 0;
  for (int i = 0; i < mx[0]; i++) {
	  for (int j = 0; j < mx[1]; j++) {
		  for (int k = 0; k < mx[2]; k++) {
              for (int l = 0; l < 3; l++) {
                  tempData[cur] = data[i][j][k][l];
                  cur++;
              }
		  }
	  }
  }

  // Write data
  H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, tempData.data());
  H5Dclose(dataset);
  H5Sclose(dataspace);
  return true;
}

template <typename T, typename>
bool Hdf5Stream::Write2DMatrix(const std::string &ArrayName, const NumMatrix<T, 2> &data,
                               const double *xb, const double *dx, bool with_opendxinfo) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

  */
  int mx[2];
  mx[0] = data.getHigh(1) - data.getLow(1) + 1;
  mx[1] = data.getHigh(0) - data.getLow(0) + 1;

  AddDatasetName(ArrayName);

  num += 1;

  const int DIM = 2;
  // Set datatype to float
  hid_t datatype = get_hdf5_data_type<T>();
  // Use little endian
  H5Tset_order(datatype, H5T_ORDER_LE);

  hsize_t DimsData[DIM];
  for (int q = 0; q < DIM; ++q) {
    DimsData[q] = mx[q];
  }
  // Make dataspace:
  hid_t dataspace = H5Screate_simple(DIM, DimsData, NULL);

  // Create dataset
  hid_t dataset = H5Dcreate2(group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);

  // Writing attributes (if necessary)
  if (with_opendxinfo) {
    // Supplying additional attributes for opendx input
    hsize_t DimsAttr = 2;
    hid_t attrspace = H5Screate_simple(1, &DimsAttr, NULL);
    hid_t datatypeAttr = get_hdf5_data_type<double>();

    double Origin[DIM];
    double Delta[DIM];
    for (int q = 0; q < DIM; ++q) {
      Origin[q] = xb[q];
      Delta[q] = dx[q];
    }

    hid_t origin = H5Acreate2(dataset, "origin", datatypeAttr, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(origin, datatypeAttr, &Origin);
    H5Aclose(origin);

    hid_t delta = H5Acreate2(dataset, "delta", datatypeAttr, attrspace, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(delta, datatypeAttr, &Delta);
    H5Aclose(delta);

    H5Sclose(attrspace);
  }

  // Write data
  H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  return true;
}

template <typename T, typename>
bool Hdf5Stream::Write3DMatrixSwap(const std::string &ArrayName, const NumMatrix<T, 3> &data,
                                   const double *xb, const double *dx) {
  /* Routine to write NumMatix data in correct ordering to hdf5 file

     Remarks: Only use if sufficient memory available - dummy
     NumMatrix has to be created
  */
  int hi[3];
  hi[0] = data.getHigh(2);
  hi[1] = data.getHigh(1);
  hi[2] = data.getHigh(0);
  int lo[3];
  lo[0] = data.getLow(2);
  lo[1] = data.getLow(1);
  lo[2] = data.getLow(0);

  NumMatrix<double, 3> outdata(Index::set(lo[0], lo[1], lo[2]), Index::set(hi[0], hi[1], hi[2]));

  // Swapping data dimensions
  for (int i = lo[0]; i <= hi[0]; ++i) {
    for (int j = lo[1]; j <= hi[1]; ++j) {
      for (int k = lo[2]; k <= hi[2]; ++k) {
        outdata(i, j, k) = data(k, j, i);
      }
    }
  }

  int mx[3];
  mx[2] = hi[0] - lo[0] + 1;
  mx[1] = hi[1] - lo[1] + 1;
  mx[0] = hi[2] - lo[2] + 1;

  num += 1;
  if (num > NumEntries) {
    return false;
  }
  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);
  const int DIM = 3;
  hsize_t DimsData[DIM];
  for (int q = 0; q < DIM; ++q) {
    DimsData[q] = mx[q];
  }
  hid_t dataspace = H5Screate_simple(DIM, DimsData, NULL);

  // Supplying additional attributes for opendx input
  // Datatype: double, little endian of size 1
  hid_t datatypefloat = get_hdf5_data_type<double>();
  H5Tset_order(datatypefloat, H5T_ORDER_LE);

  // Create dataspace for attribute
  hsize_t DimsAttr = 3;
  hid_t attrspace = H5Screate_simple(1, &DimsAttr, NULL);

  double Origin[3];
  double Delta[3];
  for (int q = 0; q < 3; ++q) {
    Origin[q] = xb[q];
    Delta[q] = dx[q];
  }

  // Create dataset
  hid_t dataset = H5Dcreate2(group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
  // Create attributes
  hid_t origin = H5Acreate2(dataset, "origin", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  hid_t delta = H5Acreate2(dataset, "delta", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);
  // Write attributes
  H5Awrite(origin, datatypefloat, &Origin);
  H5Awrite(delta, datatypefloat, &Delta);

  // Close attributes
  H5Aclose(origin);
  H5Aclose(delta);
  H5Sclose(attrspace);

  // Write data
  H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, outdata);
  H5Dclose(dataset);
  H5Sclose(dataspace);

  return true;
}

// End Hdf5Stream Templates
// -----------------------------------------------------------

// ---------------------------------------------------------
// Begin Hdf5iStream Templates

template <typename T, typename>
int Hdf5iStream::ReadNumArray(const std::string &ArrayName, NumArray<T> &data) {
  return ReadNumArray(ArrayName, data, group);
}

template <typename T, typename>
int Hdf5iStream::ReadNumArray(const std::string &ArrayName, NumArray<T> &data, hid_t my_group) {
  //! Read 1D NumArray from hdf5 file in standard group
  hid_t dataset = H5Dopen2(my_group, ArrayName.c_str(), H5P_DEFAULT);

  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset);

  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != 1) {
    std::cerr << " Input data must be one dimensional but is: " << std::endl;
    std::cerr << dimhdf << std::endl;
    exit(-2);
  }

  // Get size of dataset
  hsize_t dims_out;
  int ndims = H5Sget_simple_extent_dims(dataspace, &dims_out, NULL);

  data.resize(dims_out);

  // Read data
  int return_value = H5Dread(dataset, get_hdf5_data_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  // Close the dataspace
  H5Sclose(dataspace);
  // Close dataset
  H5Dclose(dataset);

  return return_value;
}

template <typename T, typename>
bool Hdf5iStream::ReadGlobalAttr(const std::string &AttrName, T &AttrData) {
  return ReadGlobalAttr(group, AttrName, AttrData);
}

template <typename T, typename>
bool Hdf5iStream::ReadGlobalAttr(hid_t my_group, const std::string &AttrName, T &AttrData) {
  // Open Attribute
  hid_t info = H5Aopen(my_group, AttrName.c_str(), H5P_DEFAULT);
  // Read Attribute Data
  return_val = H5Aread(info, get_hdf5_data_type<T>(), &AttrData);
  H5Aclose(info);
  return true;
}

template <typename T, typename>
int Hdf5iStream::ReadGlobalAttr(std::string GroupName, const std::string &AttrName, T &AttrData) {
  // Open group
  hid_t my_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  ReadGlobalAttr(my_group, AttrName, AttrData);
  H5Gclose(my_group);
  return return_val;
}

template <typename T, typename>
bool Hdf5iStream::Read1DMatrix(const std::string &ArrayName, NumMatrix<T, 1> &data) {
  return Read1DMatrix("Data", ArrayName, data);
}

template <typename T, typename>
bool Hdf5iStream::Read1DMatrix(std::string GroupName, const std::string &ArrayName,
                               NumMatrix<T, 1> &data) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

   */
  hid_t loc_group;
  if (GroupName != "Data") {
    loc_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  } else {
    loc_group = group;
  }

  const int DIM = 1;
  int lbound[1], ubound[1];

  hid_t dataset = H5Dopen2(loc_group, ArrayName.c_str(), H5P_DEFAULT);

  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset);

  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }
  hsize_t dims_out[DIM];
  int ndims = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

  if (ndims != DIM) {
    cerr << " Wrong number of dimensions " << ndims << " - " << DIM << endl;
    exit(-22);
  }

  lbound[0] = 0;
  ubound[0] = int(dims_out[0]) - 1;

  data.resize(lbound, ubound);

  H5Dread(dataset, get_hdf5_data_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  // Close the dataset:
  H5Dclose(dataset);
  H5Sclose(dataspace);

  if (GroupName != "Data") {
    // Close group
    H5Gclose(loc_group);
  }

  return true;
}

template <typename T, typename>
bool Hdf5iStream::Read2DMatrix(const std::string &ArrayName, NumMatrix<T, 2> &data) {
  const int DIM = 2;
  int lbound[2], ubound[2];

  hid_t dataset = H5Dopen2(group, ArrayName.c_str(), H5P_DEFAULT);

  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset);

  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }
  hsize_t dims_out[DIM];
  int ndims = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

  if (ndims != DIM) {
    cerr << " Wrong number of dimensions " << ndims << " - " << DIM << endl;
    exit(-22);
  }

  for (int i = 0; i < 2; i++) {
    lbound[i] = 0;
    ubound[i] = int(dims_out[DIM - (i + 1)]) - 1;
  }
  data.resize(lbound, ubound);
  // cout << ubound[0] << "	" << ubound[1] << endl;

  H5Dread(dataset, get_hdf5_data_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  // Close the dataset:
  H5Dclose(dataset);
  return true;
}

template <typename T, typename>
bool Hdf5iStream::Read3DMatrix(const std::string &ArrayName, NumMatrix<T, 3> &data,
                               bool load_size) {
  Read3DMatrix("Data", ArrayName, data, load_size);
  return true;
}

template <typename T, typename>
int Hdf5iStream::Read3DMatrix(std::string GroupName, const std::string &ArrayName,
                              NumMatrix<T, 3> &data, bool load_size) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

  */
  hid_t loc_group;
  if (GroupName != "Data") {
    loc_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  } else {
    loc_group = group;
  }

  const int DIM = 3;
  int mx[3];
  if (!load_size) {
    mx[0] = data.getHigh(2) - data.getLow(2) + 1;
    mx[1] = data.getHigh(1) - data.getLow(1) + 1;
    mx[2] = data.getHigh(0) - data.getLow(0) + 1;
  }

  hid_t dataset = H5Dopen2(loc_group, ArrayName.c_str(), H5P_DEFAULT);

  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset); /* dataspace handle */

  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }
  hsize_t dims_out[DIM];
  // Get dims
  int ndims = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
  if (ndims != DIM) {
    cerr << " Wrong number of dimensions " << ndims << " - " << DIM << endl;
    exit(-22);
  }

  cout << " Trying to read " << ArrayName << " " << load_size << endl;

  if (load_size) {
    int lbound[3], ubound[3];
    for (int i = 0; i < 3; i++) {
      lbound[i] = 0;
      ubound[i] = int(dims_out[DIM - (i + 1)]) - 1;
    }
    data.resize(lbound, ubound);
  } else {
    for (int i = 0; i < DIM; ++i) {
      if ((int(dims_out[i])) != mx[i]) {
        cerr << " Wrong size of dimension " << i << ":" << endl;
        cerr << int(dims_out[i]) << " " << mx[i] << endl;
        exit(-33);
      }
    }
  }

  int q_index = -1;
  if (H5Aexists(dataset, "q_index")) {
    hid_t attr = H5Aopen(dataset, "q_index", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &q_index);
    H5Aclose(attr);
  }

  H5Dread(dataset, get_hdf5_data_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  // close dataset
  H5Dclose(dataset);

  if (GroupName != "Data") {
    // Close group
    H5Gclose(loc_group);
  }

  return q_index;
}

template <typename T, typename>
bool Hdf5iStream::Read2DFrom3DMatrix(const std::string &ArrayName, NumMatrix<T, 2> &data, int dir,
                                     int posPerp) {
  return Read2DFrom3DMatrix("Data", ArrayName, data, dir, posPerp);
}

template <typename T, typename>
bool Hdf5iStream::Read2DFrom3DMatrix(std::string GroupName, const std::string &ArrayName,
                                     NumMatrix<T, 2> &data, int dir, int posPerp) {
  hid_t loc_group;
  if (GroupName != "Data") {
    loc_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  } else {
    loc_group = group;
  }

  // Open the dataset:
  hid_t dataset = H5Dopen2(loc_group, ArrayName.c_str(), H5P_DEFAULT);

  // Get dimensions of dataset:
  const int DIM = 3;

  hid_t dataspace = H5Dget_space(dataset);
  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }
  hsize_t dims_out[DIM];
  int ndims = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

  // Beware - dimensions are swapped
  hsize_t mx[2];
  if (dir == 0) {  // x,y
    mx[0] = dims_out[1];
    mx[1] = dims_out[2];
  } else if (dir == 1) {  // x,z
    mx[0] = dims_out[0];
    mx[1] = dims_out[2];
  } else if (dir == 2) {  // y,z
    mx[0] = dims_out[0];
    mx[1] = dims_out[1];
  } else {
    cerr << " Error no such direction " << endl;
    exit(-234);
  }

  cerr << " Making memspace " << endl;
  hid_t memspace = H5Screate_simple(2, mx, NULL);
  cerr << " done " << endl;

  hsize_t count[3], offset[3];

  int lbound[2], ubound[2];
  lbound[0] = 0;
  lbound[1] = 0;

  if (dir == 0) {  // x,y plane
    // count[0] = mx[0];
    // count[1] = mx[1];
    // count[2] = 1;
    // offset[0] = 0;
    // offset[1] = 0;
    // offset[2] = posPerp;
    count[0] = 1;
    count[1] = dims_out[1];
    count[2] = dims_out[2];
    offset[0] = posPerp;
    offset[1] = 0;
    offset[2] = 0;
    ubound[0] = dims_out[2] - 1;
    ubound[1] = dims_out[1] - 1;
  } else if (dir == 1) {  // x,z plane
    // count[0] = mx[0];
    // count[1] = 1;
    // count[2] = mx[1];
    // offset[0] = 0;
    // offset[1] = posPerp;
    // offset[2] = 0;
    count[0] = dims_out[0];
    count[1] = 1;
    count[2] = dims_out[2];
    offset[0] = 0;
    offset[1] = posPerp;
    offset[2] = 0;
    ubound[0] = dims_out[2] - 1;
    ubound[1] = dims_out[0] - 1;
  } else {  // y,z plane
    // count[0] = 1;
    // count[1] = mx[0];
    // count[2] = mx[1];
    // offset[0] = posPerp;
    // offset[1] = 0;
    // offset[2] = 0;
    count[0] = dims_out[0];
    count[1] = dims_out[1];
    count[2] = 1;
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = posPerp;
    ubound[0] = dims_out[1] - 1;
    ubound[1] = dims_out[0] - 1;
  }

  cerr << " resize " << endl;
  data.resize(lbound, ubound);
  cerr << " resize done " << endl;

  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
  cerr << " selc hyper done " << endl;
  H5Dread(dataset, get_hdf5_data_type<T>(), memspace, dataspace, H5P_DEFAULT, data);
  cerr << " read done " << endl;

  // Close the dataset:
  H5Dclose(dataset);
  H5Sclose(memspace);
  if (GroupName != "Data") {
    // Close group
    H5Gclose(loc_group);
  }
  return true;
}

template <typename T, typename>
T Hdf5iStream::ReadPointFrom3DMatrix(const std::string &ArrayName, int ix, int iy, int iz) {
  // Open the dataset:
  hid_t dataset = H5Dopen2(group, ArrayName.c_str(), H5P_DEFAULT);

  // Get dimensions of dataset:
  const int DIM = 3;

  hid_t dataspace = H5Dget_space(dataset);
  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }

  // Select just a single point
  hsize_t dim_point[] = {1};

  hid_t memspace = H5Screate_simple(1, dim_point, NULL);

  const int numPoints = 1;
  // Set coordinates of a point
  // hssize_t coord[numPoints][DIM];
  // coord[0][0] = iz;
  // coord[0][1] = iy;
  // coord[0][2] = ix;
  hsize_t coord[DIM];  // = {iz, iy, ix};
  coord[0] = iz;
  coord[1] = iy;
  coord[2] = ix;

  // H5Sselect_elements(dataspace, H5S_SELECT_SET, numPoints,
  //                    (const hssize_t **)coord);
  H5Sselect_elements(dataspace, H5S_SELECT_SET, 1, coord);

  float data[numPoints];
  H5Dread(dataset, get_hdf5_data_type<T>(), memspace, dataspace, H5P_DEFAULT, data);

  cout << " Data: " << data[0] << endl;

  // Close the dataset:
  H5Dclose(dataset);
  H5Sclose(memspace);

  return data[0];
}
template <typename T, typename>
bool Hdf5iStream::Read3DMatrix(const std::string &ArrayName, NumMatrix<T, 3> &data, double *xb,
                               double *dx) {
  return Read3DMatrix("Data", ArrayName, data, xb, dx);
}
template <typename T, typename>
bool Hdf5iStream::Read3DMatrix(std::string GroupName, const std::string &ArrayName,
                               NumMatrix<T, 3> &data, double *xb, double *dx) {
  /*
     Routine to read 3D matrix from h5-file. Parameters are:
     GroupName -> Name of the data-group
     ArrayName -> Name of dataset w/o group name
     data      -> NumMatrix-Array to hold data
     xb        -> Array to hold values of lower bound positions
     dx        -> Array to hold grid size
  */

  const int DIM = 3;
  int lbound[3], ubound[3];
  float dummy[3];

  hid_t loc_group;
  if (GroupName != "Data") {
    loc_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  } else {
    loc_group = group;
  }

  hid_t dataset = H5Dopen2(loc_group, ArrayName.c_str(), H5P_DEFAULT);
  hid_t Origin = H5Aopen(dataset, "origin", H5P_DEFAULT);
  H5Aread(Origin, H5T_NATIVE_DOUBLE, dummy);
  xb[0] = double(dummy[0]);
  xb[1] = double(dummy[1]);
  xb[2] = double(dummy[2]);

  hid_t Delta = H5Aopen(dataset, "delta", H5P_DEFAULT);
  H5Aread(Delta, H5T_NATIVE_DOUBLE, dummy);
  dx[0] = double(dummy[0]);
  dx[1] = double(dummy[1]);
  dx[2] = double(dummy[2]);

  // Get dataspace handler
  hid_t dataspace = H5Dget_space(dataset);

  int dimhdf = H5Sget_simple_extent_ndims(dataspace);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }
  hsize_t dims_out[DIM];
  H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

  for (int i = 0; i < 3; i++) {
    lbound[i] = 0;
    ubound[i] = int(dims_out[DIM - (i + 1)]) - 1;
  }
  data.resize(lbound, ubound);

  H5Dread(dataset, get_hdf5_data_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  // Close attribute ids
  H5Aclose(Origin);
  H5Aclose(Delta);

  // Close the dataset:
  H5Dclose(dataset);

  if (GroupName != "Data") {
    // Close group
    H5Gclose(loc_group);
  }

  return true;
}

template <typename T, typename>
T Hdf5iStream::ReadDatasetAttr(std::string DataSetName, const std::string &AttrName) {
  // Open the dataset
  hid_t dataset = H5Dopen2(group, DataSetName.c_str(), H5P_DEFAULT);

  hid_t attr = H5Aopen(dataset, AttrName.c_str(), H5P_DEFAULT);
  T value;
  H5Aread(attr, get_hdf5_data_type<T>(), &value);

  H5Aclose(attr);
  H5Dclose(dataset);

  return value;
}

template <typename T, typename>
int Hdf5iStream::Read3DMatrix_parallel(const std::string &ArrayName, const NumArray<int> &mx_local,
                                       const NumArray<int> &rank_shift, NumMatrix<T, 3> &data) {
  return Read3DMatrix_parallel("Data", ArrayName, mx_local, rank_shift, data);
}

template <typename T, typename>
int Hdf5iStream::Read3DMatrix_parallel(std::string GroupName, const std::string &ArrayName,
                                       const NumArray<int> &mx_local,
                                       const NumArray<int> &rank_shift, NumMatrix<T, 3> &data) {
  /* Routine to write NumMatrix data in wrong ordering to hdf5 file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

  */
  hid_t loc_group;
  if (GroupName != "Data") {
    loc_group = H5Gopen2(hdf5file, GroupName.c_str(), H5P_DEFAULT);
  } else {
    loc_group = group;
  }

  // Determine rim of data
  int rim = 0;
  if (data.getLow(0) < 0) {
    rim = -data.getLow(0);
  }

  const int DIM = 3;
  int mx[3];
  mx[0] = data.getHigh(2) - data.getLow(2) - 2 * rim;
  mx[1] = data.getHigh(1) - data.getLow(1) - 2 * rim;
  mx[2] = data.getHigh(0) - data.getLow(0) - 2 * rim;

  hid_t dataset_id = H5Dopen2(loc_group, ArrayName.c_str(), H5P_DEFAULT);

  // Get dataspace handler
  hid_t dataspace_id = H5Dget_space(dataset_id); /* dataspace handle */

  cout << " Trying to read " << ArrayName << endl;

  int dimhdf = H5Sget_simple_extent_ndims(dataspace_id);
  if (dimhdf != DIM) {
    cerr << " Wrong dimensionality of input data: " << endl;
    cerr << dimhdf << " " << DIM << endl;
    exit(-2);
  }
  hsize_t dims_out[DIM];
  // Get dims
  int ndims = H5Sget_simple_extent_dims(dataspace_id, dims_out, NULL);
  if (ndims != DIM) {
    cerr << " Wrong number of dimensions " << ndims << " - " << DIM << endl;
    exit(-22);
  }
  // Here, we only have to check that sufficient data is available
  for (int i = 0; i < DIM; ++i) {
    if (mx[i] > (int(dims_out[i]))) {
      cerr << " Wrong size of dimension " << i << ":" << endl;
      cerr << int(dims_out[i]) << " " << mx[i] << endl;
    }
  }

  // Check if mx_local is compatible with mx
  for (int dir = 0; dir < DIM; ++dir) {
    if (mx[DIM - dir - 1] != mx_local[dir]) {
      cerr << " Wrong size of dimension " << dir << ":" << endl;
      cerr << int(dims_out[dir]) << " " << mx[dir] << endl;
      exit(3);
    }
  }

  // Compute shift of local dataset w.r.t. data in file
  hsize_t offset[DIM];
  hsize_t count_local[DIM];
  for (int i_dir = 0; i_dir < DIM; ++i_dir) {
    offset[i_dir] = rank_shift(DIM - i_dir - 1);
    count_local[i_dir] = mx_local(DIM - i_dir - 1) + 1 + 2 * rim;
  }

  // Create local dataspace
  hid_t dataspaceLocal = H5Screate_simple(dimhdf, count_local, NULL);

  // Select local data as hyperslab in dataset
  H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count_local, NULL);

  int q_index = -1;
  if (H5Aexists(dataset_id, "q_index")) {
    hid_t attr = H5Aopen(dataset_id, "q_index", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &q_index);
    H5Aclose(attr);
  }
  // Read data to local array
  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  hid_t return_val =
      H5Dread(dataset_id, get_hdf5_data_type<T>(), dataspaceLocal, dataspace_id, plist_id, data);

  // Close property list
  H5Pclose(plist_id);

  // Close local dataspace
  H5Sclose(dataspaceLocal);

  // close dataset
  H5Dclose(dataset_id);

  if (GroupName != "Data") {
    // Close group
    H5Gclose(loc_group);
  }

  return q_index;
}

#if (HDF_PARALLEL_IO == CRONOS_ON)
template <typename T, typename>
bool Hdf5iStream::Read3DMatrix_withMPI_IO(const std::string &ArrayName, NumMatrix<T, 3> &data,
                                          NumArray<int> &mx_local, NumArray<int> &rank_shift) {
  return Read3DMatrix_withMPI_IO(ArrayName, data, mx_local, rank_shift, group);
}

template <typename T, typename>
bool Hdf5iStream::Read3DMatrix_withMPI_IO(const std::string &ArrayName, NumMatrix<T, 3> &data,
                                          NumArray<int> &mx_local, NumArray<int> &rank_shift,
                                          std::string GroupName) {
  hid_t my_group = OpenGroup(GroupName);
  bool ret = Read3DMatrix_withMPI_IO(ArrayName, data, mx_local, rank_shift, my_group);
  CloseGroup(my_group);
  return ret;
}

template <typename T, typename>
bool Hdf5iStream::Read3DMatrix_withMPI_IO(const std::string &ArrayName, NumMatrix<T, 3> &data,
                                          NumArray<int> &mx_local, NumArray<int> &rank_shift,
                                          hid_t my_group) {
  /* Routine to read NumMatrix data in wrong ordering to hdf5 file. MPI
     parallel form of the routine that reads all data from a single file.

     Remarks:

     On reading the hdf data the swapped dimensions have to be taken
     into account

     rim can be computed from mx_local.

  */

  int DIM = 3;
  float dummy[3];

  // Determine rim of data
  int rim = 0;
  if (data.getLow(0) < 0) {
    rim = -data.getLow(0);
  }

  int mx[3];
  mx[0] = data.getHigh(2) - data.getLow(2) + 1;
  mx[1] = data.getHigh(1) - data.getLow(2) + 1;
  mx[2] = data.getHigh(0) - data.getLow(2) + 1;

  // open the dataset
  hid_t dataset = H5Dopen2(my_group, ArrayName.c_str(), H5P_DEFAULT);

  hid_t datatype = get_hdf5_data_type<T>();
  H5Tset_order(datatype, H5T_ORDER_LE);

  hid_t dataspace = H5Dget_space(dataset);

  // Now make distinguis local from global data via hyperslabs:
  hsize_t sizeLocal[DIM];
  sizeLocal[0] = mx_local[2] + 1 + 2 * rim;
  sizeLocal[1] = mx_local[1] + 1 + 2 * rim;
  sizeLocal[2] = mx_local[0] + 1 + 2 * rim;

  hsize_t offset[DIM];
  offset[0] = rank_shift[2];
  offset[1] = rank_shift[1];
  offset[2] = rank_shift[0];

  hid_t dataspaceLocal = H5Screate_simple(DIM, sizeLocal, NULL);

  // Select local data as hyperslab in dataset
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, sizeLocal, NULL);

  // Read data
  H5Dread(dataset, datatype, dataspaceLocal, dataspace, H5P_DEFAULT, data);

  H5Sclose(dataspaceLocal);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return true;
}
#endif


// End Hdf5iStream Templates
// -----------------------------------------------------------

// ---------------------------------------------------------
// Evaluate Templates for the most common datatypes

template bool Hdf5Stream::ChangeGlobalAttr(const std::string &AttrName, bool AttrData);
template bool Hdf5Stream::ChangeGlobalAttr(const std::string &AttrName, int AttrData);
template bool Hdf5Stream::ChangeGlobalAttr(const std::string &AttrName, float AttrData);
template bool Hdf5Stream::ChangeGlobalAttr(const std::string &AttrName, double AttrData);

template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, bool AttrData);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, std::string AttrData);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, int AttrData);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, float AttrData);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, double AttrData);

//template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, bool AttrData, hid_t my_group);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, std::string AttrData,
                                        hid_t my_group);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, int AttrData, hid_t my_group);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, float AttrData,
                                        hid_t my_group);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, double AttrData,
                                        hid_t my_group);

template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const char *AttrData,
                                        int entries);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const bool *AttrData,
                                        int entries);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const std::string *AttrData,
                                        int entries);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const int *AttrData,
                                        int entries);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const float *AttrData,
                                        int entries);
template bool Hdf5Stream::AddGlobalAttr(const std::string &AttrName, const double *AttrData,
                                        int entries);

template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName, bool AttributeData);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              std::string AttributeData);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName, int AttributeData);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              float AttributeData);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              double AttributeData);

template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName, bool AttributeData,
                                              hid_t my_group);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              std::string AttributeData, hid_t my_group);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName, int AttributeData,
                                              hid_t my_group);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName, float AttributeData,
                                              hid_t my_group);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              double AttributeData, hid_t my_group);

template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              const double *AttributeData, const hsize_t entries);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              const float *AttributeData, const hsize_t entries);
template bool Hdf5Stream::AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              const int *AttributeData, const hsize_t entries);

template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName, const NumMatrix<int, 1> &data,
                                        double Origin, double Delta, int numin);
template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName,
                                        const NumMatrix<float, 1> &data, double Origin,
                                        double Delta, int numin);
template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName,
                                        const NumMatrix<double, 1> &data, double Origin,
                                        double Delta, int numin);

template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName,
                                        const NumMatrix<int, 1> &data);
template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName,
                                        const NumMatrix<float, 1> &data);
template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName,
                                        const NumMatrix<double, 1> &data);

template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName, const NumMatrix<int, 1> &data,
                                        hid_t my_group);
template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName,
                                        const NumMatrix<float, 1> &data, hid_t my_group);
template bool Hdf5Stream::Write1DMatrix(const std::string &ArrayName,
                                        const NumMatrix<double, 1> &data, hid_t my_group);

template bool Hdf5Stream::WriteNumArray(const std::string &ArrayName, const NumArray<int> &data);
template bool Hdf5Stream::WriteNumArray(const std::string &ArrayName, const NumArray<float> &data);
template bool Hdf5Stream::WriteNumArray(const std::string &ArrayName, const NumArray<double> &data);

template bool Hdf5Stream::WriteNumArray(const std::string &ArrayName, const NumArray<int> &data,
                                        hid_t my_group);
template bool Hdf5Stream::WriteNumArray(const std::string &ArrayName, const NumArray<float> &data,
                                        hid_t my_group);
template bool Hdf5Stream::WriteNumArray(const std::string &ArrayName, const NumArray<double> &data,
                                        hid_t my_group);

template bool Hdf5Stream::Write3DMatrix(const std::string &ArrayName,
                                        const NumMatrix<float, 3> &data, const double *xb,
                                        const double *dx, hid_t my_group, int q_index,
                                        bool with_opendxinfo);

template bool Hdf5Stream::Write3DMatrix(const std::string &ArrayName,
                                        const NumMatrix<double, 3> &data, const double *xb,
                                        const double *dx, hid_t my_group, int q_index,
                                        bool with_opendxinfo);

template bool Hdf5Stream::Write3DMatrix(const std::string &ArrayName,
                                        const NumMatrix<float, 3> &data, const double *xb,
                                        const double *dx);
template bool Hdf5Stream::Write3DMatrix(const std::string &ArrayName,
                                        const NumMatrix<double, 3> &data, const double *xb,
                                        const double *dx);

template bool Hdf5Stream::Write3DMatrix(const std::string &ArrayName,
                                        const NumMatrix<float, 3> &data);

template bool Hdf5Stream::WriteNDArray(const std::string &ArrayName, const float *data,
                                       const int mx[], int dim);

template bool Hdf5Stream::Write3DVecMatrix(const std::string &ArrayName,
                                           const NumMatrix<float, 3> &data_x,
                                           const NumMatrix<float, 3> &data_y,
                                           const NumMatrix<float, 3> &data_z);

template bool Hdf5Stream::Write2DMatrix(const std::string &ArrayName,
                                        const NumMatrix<float, 2> &data, const double *xb,
                                        const double *dx, bool with_opendxinfo /*= true*/);
template bool Hdf5Stream::Write2DMatrix(const std::string &ArrayName,
                                        const NumMatrix<double, 2> &data, const double *xb,
                                        const double *dx, bool with_opendxinfo /*= true*/);

template bool Hdf5Stream::Write3DMatrixSwap(const std::string &ArrayName,
                                            const NumMatrix<double, 3> &data, const double *xb,
                                            const double *dx);


#if (HDF_PARALLEL_IO == CRONOS_ON)
template bool Hdf5Stream::WriteNumArray_withMPI_IO(const std::string &ArrayName,
                                                   const NumArray<float> &data, bool with_data);
template bool Hdf5Stream::WriteNumArray_withMPI_IO(const std::string &ArrayName,
                                                   const NumArray<double> &data, bool with_data);
template bool Hdf5Stream::WriteNumArray_withMPI_IO2(const std::string &ArrayName,
                                                    const NumArray<float> &data, bool with_data);
template bool Hdf5Stream::Write3DMatrix_withMPI_IO(
    const std::string &ArrayName, const NumMatrix<float, 3> &data, const NumArray<int> &mx_global,
    const NumArray<int> &mx_local, const NumArray<int> &rank_shift, const NumArray<float> &xb,
    const NumArray<float> &dx);
template bool Hdf5Stream::Write3DMatrix_withMPI_IO(
    const std::string &ArrayName, const NumMatrix<float, 3> &data, const NumArray<int> &mx_global,
    const NumArray<int> &mx_local, const NumArray<int> &rank_shift,
    // NumArray<int> &rank_pos,
    const NumArray<float> &xb, const NumArray<float> &dx, hid_t my_group, int q_index,
    bool with_opendxinfo);
template bool Hdf5Stream::Write3DMatrix_withMPI_IO(
    const std::string &ArrayName, const NumMatrix<double, 3> &data, const NumArray<int> &mx_global,
    const NumArray<int> &mx_local, const NumArray<int> &rank_shift,
    // NumArray<int> &rank_pos,
    const NumArray<float> &xb, const NumArray<float> &dx, hid_t my_group, int q_index,
    bool with_opendxinfo);

template bool Hdf5Stream::Write2DMatrix_withMPI_IO(
    const std::string &ArrayName, const NumMatrix<double, 2> &data, const NumArray<int> &mx_global,
    const NumArray<int> &mx_local, const NumArray<int> &rank_shift, const NumArray<float> &xb,
    const NumArray<float> &dx, hid_t my_group, bool);

template bool Hdf5Stream::Write2DMatrix_withMPI_IO(
    const std::string &ArrayName, const NumMatrix<double, 2> &data, const NumArray<int> &mx_global,
    const NumArray<int> &mx_local, const NumArray<int> &rank_shift, const NumArray<float> &xb,
    const NumArray<float> &dx);
#endif


// ----------------------------------------------------------------------------------------

template int Hdf5iStream::ReadNumArray(const std::string &ArrayName, NumArray<int> &data);
template int Hdf5iStream::ReadNumArray(const std::string &ArrayName, NumArray<float> &data);
template int Hdf5iStream::ReadNumArray(const std::string &ArrayName, NumArray<double> &data);

template int Hdf5iStream::ReadNumArray(const std::string &ArrayName, NumArray<int> &data,
                                       hid_t my_group);
template int Hdf5iStream::ReadNumArray(const std::string &ArrayName, NumArray<float> &data,
                                       hid_t my_group);
template int Hdf5iStream::ReadNumArray(const std::string &ArrayName, NumArray<double> &data,
                                       hid_t my_group);

template bool Hdf5iStream::ReadGlobalAttr(const std::string &AttrName, bool &AttrData);
template bool Hdf5iStream::ReadGlobalAttr(const std::string &AttrName, std::string &AttrData);
template bool Hdf5iStream::ReadGlobalAttr(const std::string &AttrName, int &AttrData);
template bool Hdf5iStream::ReadGlobalAttr(const std::string &AttrName, float &AttrData);
template bool Hdf5iStream::ReadGlobalAttr(const std::string &AttrName, double &AttrData);

//template bool Hdf5iStream::ReadGlobalAttr(hid_t my_group, const std::string &AttrName,
//                                          bool &AttrData);
template bool Hdf5iStream::ReadGlobalAttr(hid_t my_group, const std::string &AttrName,
                                          std::string &AttrData);
template bool Hdf5iStream::ReadGlobalAttr(hid_t my_group, const std::string &AttrName,
                                          int &AttrData);
template bool Hdf5iStream::ReadGlobalAttr(hid_t my_group, const std::string &AttrName,
                                          float &AttrData);
template bool Hdf5iStream::ReadGlobalAttr(hid_t my_group, const std::string &AttrName,
                                          double &AttrData);

template int Hdf5iStream::ReadGlobalAttr(std::string GroupName, const std::string &AttrName,
                                         bool &AttrData);
template int Hdf5iStream::ReadGlobalAttr(std::string GroupName, const std::string &AttrName,
                                         std::string &AttrData);
template int Hdf5iStream::ReadGlobalAttr(std::string GroupName, const std::string &AttrName,
                                         int &AttrData);
template int Hdf5iStream::ReadGlobalAttr(std::string GroupName, const std::string &AttrName,
                                         float &AttrData);
template int Hdf5iStream::ReadGlobalAttr(std::string GroupName, const std::string &AttrName,
                                         double &AttrData);

template bool Hdf5iStream::Read1DMatrix(const std::string &ArrayName, NumMatrix<int, 1> &data);
template bool Hdf5iStream::Read1DMatrix(const std::string &ArrayName, NumMatrix<float, 1> &data);
template bool Hdf5iStream::Read1DMatrix(const std::string &ArrayName, NumMatrix<double, 1> &data);

template bool Hdf5iStream::Read1DMatrix(std::string GroupName, const std::string &ArrayName,
                                        NumMatrix<int, 1> &data);
template bool Hdf5iStream::Read1DMatrix(std::string GroupName, const std::string &ArrayName,
                                        NumMatrix<float, 1> &data);
template bool Hdf5iStream::Read1DMatrix(std::string GroupName, const std::string &ArrayName,
                                        NumMatrix<double, 1> &data);

template bool Hdf5iStream::Read3DMatrix(const std::string &ArrayName, NumMatrix<float, 3> &data,
                                        bool load_size);
template bool Hdf5iStream::Read3DMatrix(const std::string &ArrayName, NumMatrix<double, 3> &data,
                                        bool load_size);

template int Hdf5iStream::Read3DMatrix(std::string GroupName, const std::string &ArrayName,
                                       NumMatrix<float, 3> &data, bool load_size);
template int Hdf5iStream::Read3DMatrix(std::string GroupName, const std::string &ArrayName,
                                       NumMatrix<double, 3> &data, bool load_size);

template bool Hdf5iStream::Read3DMatrix(const std::string &ArrayName, NumMatrix<float, 3> &data,
                                        double *xb, double *dx);
template bool Hdf5iStream::Read3DMatrix(const std::string &ArrayName, NumMatrix<double, 3> &data,
                                        double *xb, double *dx);

template bool Hdf5iStream::Read3DMatrix(std::string GroupName, const std::string &ArrayName,
                                        NumMatrix<float, 3> &data, double *xb, double *dx);
template bool Hdf5iStream::Read3DMatrix(std::string GroupName, const std::string &ArrayName,
                                        NumMatrix<double, 3> &data, double *xb, double *dx);

template bool Hdf5iStream::Read2DMatrix(const std::string &ArrayName, NumMatrix<int, 2> &);
template bool Hdf5iStream::Read2DMatrix(const std::string &ArrayName, NumMatrix<float, 2> &);
template bool Hdf5iStream::Read2DMatrix(const std::string &ArrayName, NumMatrix<double, 2> &);

template bool Hdf5iStream::Read2DFrom3DMatrix(std::string GroupName, const std::string &ArrayName,
                                              NumMatrix<int, 2> &data, int dir, int posPerp);
template bool Hdf5iStream::Read2DFrom3DMatrix(std::string GroupName, const std::string &ArrayName,
                                              NumMatrix<float, 2> &data, int dir, int posPerp);
template bool Hdf5iStream::Read2DFrom3DMatrix(std::string GroupName, const std::string &ArrayName,
                                              NumMatrix<double, 2> &data, int dir, int posPerp);
template bool Hdf5iStream::Read2DFrom3DMatrix(const std::string &ArrayName, NumMatrix<int, 2> &data,
                                              int dir, int posPerp);
template bool Hdf5iStream::Read2DFrom3DMatrix(const std::string &ArrayName,
                                              NumMatrix<float, 2> &data, int dir, int posPerp);
template bool Hdf5iStream::Read2DFrom3DMatrix(const std::string &ArrayName,
                                              NumMatrix<double, 2> &data, int dir, int posPerp);

template float Hdf5iStream::ReadPointFrom3DMatrix(const std::string &ArrayName, int ix, int iy,
                                                  int iz);

template float Hdf5iStream::ReadDatasetAttr(std::string DataSetName, const std::string &AttrName);

template int Hdf5iStream::Read3DMatrix_parallel(const std::string &ArrayName,
                                                const NumArray<int> &mx_local,
                                                const NumArray<int> &rank_shift,
                                                NumMatrix<int, 3> &);
template int Hdf5iStream::Read3DMatrix_parallel(const std::string &ArrayName,
                                                const NumArray<int> &mx_local,
                                                const NumArray<int> &rank_shift,
                                                NumMatrix<float, 3> &);
template int Hdf5iStream::Read3DMatrix_parallel(const std::string &ArrayName,
                                                const NumArray<int> &mx_local,
                                                const NumArray<int> &rank_shift,
                                                NumMatrix<double, 3> &);

template int Hdf5iStream::Read3DMatrix_parallel(std::string GroupName, const std::string &ArrayName,
                                                const NumArray<int> &mx_local,
                                                const NumArray<int> &rank_shift,
                                                NumMatrix<int, 3> &);
template int Hdf5iStream::Read3DMatrix_parallel(std::string GroupName, const std::string &ArrayName,
                                                const NumArray<int> &mx_local,
                                                const NumArray<int> &rank_shift,
                                                NumMatrix<float, 3> &);
template int Hdf5iStream::Read3DMatrix_parallel(std::string GroupName, const std::string &ArrayName,
                                                const NumArray<int> &mx_local,
                                                const NumArray<int> &rank_shift,
                                                NumMatrix<double, 3> &);

#if (HDF_PARALLEL_IO == CRONOS_ON)
template bool Hdf5iStream::Read3DMatrix_withMPI_IO(const std::string &ArrayName,
                                                   NumMatrix<float, 3> &data,
                                                   NumArray<int> &mx_local,
                                                   NumArray<int> &rank_shift);
template bool Hdf5iStream::Read3DMatrix_withMPI_IO(const std::string &ArrayName,
                                                   NumMatrix<double, 3> &data,
                                                   NumArray<int> &mx_local,
                                                   NumArray<int> &rank_shift);

template bool Hdf5iStream::Read3DMatrix_withMPI_IO(const std::string &ArrayName,
                                                   NumMatrix<float, 3> &data,
                                                   NumArray<int> &mx_local,
                                                   NumArray<int> &rank_shift,
                                                   std::string GroupName);
template bool Hdf5iStream::Read3DMatrix_withMPI_IO(const std::string &ArrayName,
                                                   NumMatrix<double, 3> &data,
                                                   NumArray<int> &mx_local,
                                                   NumArray<int> &rank_shift,
                                                   std::string GroupName);
#endif

