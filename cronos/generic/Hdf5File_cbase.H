#ifndef CRONOS_HDF5FILE_H
#define CRONOS_HDF5FILE_H = 1

#include <stdlib.h>
#include <fstream>
#include <stdexcept>
#include <type_traits>
#include "cronos_switches.H"
#include "hdf5.h"
#include "matrix.H"

// template<typename T> hid_t get_hdf5_data_type();

// Some definitions for convenience
template <typename T>
using Invoke = typename T::type;
template <typename Condition>
using negation = std::integral_constant<bool, !bool(Condition::value)>;
template <typename Condition>
using EnableIf = Invoke<std::enable_if<Condition::value>>;
template <typename Condition>
using DisableIf = EnableIf<::negation<Condition>>;

// Do not deduce from templates for pointers
#define template_noPointer template <typename T, typename = DisableIf<std::is_pointer<T>>>

// Here the user can choose if a parallel IO to a single hdf5 file is
// desired:
//#define HDF_PARALLEL_IO CRONOS_ON

class Hdf5Stream {
 public:
  Hdf5Stream(std::string, int, int rank = 0, bool use_MPI_IO = false
#ifdef parallel
#if (HDF_PARALLEL_IO == CRONOS_ON)
             ,
             MPI_Comm = MPI_COMM_WORLD
#endif
#endif
  );
  ~Hdf5Stream();
  bool close();
  void Reopen();

  template_noPointer bool ChangeGlobalAttr(const std::string &AttrName, T AttrData);
  template_noPointer bool ChangeGlobalAttr(const std::string &AttrName, T AttrData, hid_t my_group);
  bool AddGlobalAttr(const std::string &AttrName, const char *AttrData);
  template_noPointer bool AddGlobalAttr(const std::string &AttrName, T AttrData);
  template_noPointer bool AddGlobalAttr(const std::string &AttrName, T AttrData, hid_t my_group);
  template_noPointer bool AddGlobalAttr(const std::string &AttrName, const T *AttrData,
                                        int entries);
  template_noPointer bool AddGlobalAttr(const std::string &AttrName, const T *AttrData,
                                        hid_t my_group, int entries);
  template_noPointer bool AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName, T AttributeData);
  template_noPointer bool AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName, T AttributeData,
                                              hid_t my_group);
  template_noPointer bool AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              const T *AttributeData, const hsize_t entries);
  template_noPointer bool AddAttributeToArray(const std::string &ArrayName,
                                              const std::string &AttributeName,
                                              const T *AttributeData, const hsize_t entries,
                                              hid_t my_group);
  template_noPointer bool AddAttrToArrSingle(const std::string &ArrayName, hid_t my_group,
                                             const std::string &AttributeName, T AttributeData);
  template_noPointer bool Write1DMatrix(const std::string &ArrayName, const NumMatrix<T, 1> &data,
                                        double Origin, double Delta, int numin);
  template_noPointer bool Write1DMatrix(const std::string &ArrayName, const NumMatrix<T, 1> &data);
  template_noPointer bool Write1DMatrix(const std::string &ArrayName, const NumMatrix<T, 1> &data,
                                        hid_t my_group);
  template_noPointer bool WriteNumArray(const std::string &ArrayName, const NumArray<T> &data);
  template_noPointer bool WriteNumArray(const std::string &ArrayName, const NumArray<T> &data,
                                        hid_t my_group);

  template_noPointer bool Write3DMatrix(const std::string &ArrayName, const NumMatrix<T, 3> &data,
                                        const double *xb, const double *dx);
  template_noPointer bool Write3DMatrix(const std::string &ArrayName, const NumMatrix<T, 3> &data,
                                        const double *xb, const double *dx, hid_t my_group,
                                        int q_index = -1, bool with_opendxinfo = true);

  template_noPointer bool Write3DMatrix(const std::string &ArrayName, const NumMatrix<T, 3> &data);
  template_noPointer bool WriteNDArray(const std::string &ArrayName, const T *data, const int mx[],
                                       int dim);
  template_noPointer bool Write3DVecMatrix(const std::string &ArrayName,
                                           const NumMatrix<T, 3> &data_x,
                                           const NumMatrix<T, 3> &data_y,
                                           const NumMatrix<T, 3> &data_z);
  template_noPointer bool Write2DMatrix(const std::string &ArrayName, const NumMatrix<T, 2> &data,
                                        const double *xb, const double *dx,
                                        bool with_opendxinfo = true);
  template_noPointer bool Write3DMatrixSwap(const std::string &ArrayName,
                                            const NumMatrix<T, 3> &data, const double *xb,
                                            const double *dx);

#ifdef parallel
#if (HDF_PARALLEL_IO == CRONOS_ON)
  template_noPointer bool WriteNumArray_withMPI_IO(const std::string &ArrayName,
                                                   const NumArray<T> &data, bool with_data);
  template_noPointer bool WriteNumArray_withMPI_IO2(const std::string &ArrayName,
                                                    const NumArray<T> &data, bool with_data);

  template_noPointer bool Write3DMatrix_withMPI_IO(
      const std::string &ArrayName, const NumMatrix<T, 3> &data, const NumArray<int> &mx_global,
      const NumArray<int> &mx_local, const NumArray<int> &rank_shift,
      // NumArray<int> &rank_pos,
      const NumArray<float> &xb, const NumArray<float> &dx);
  template_noPointer bool Write3DMatrix_withMPI_IO(
      const std::string &ArrayName, const NumMatrix<T, 3> &data, const NumArray<int> &mx_global,
      const NumArray<int> &mx_local, const NumArray<int> &rank_shift,
      // NumArray<int> &rank_pos,
      const NumArray<float> &xb, const NumArray<float> &dx, hid_t my_group, int q_index = -1,
      bool with_opendxinfo = true);

  template_noPointer bool Write2DMatrix_withMPI_IO(
      const std::string &ArrayName, const NumMatrix<T, 2> &data, const NumArray<int> &mx_global,
      const NumArray<int> &mx_local, const NumArray<int> &rank_shift, const NumArray<float> &xb,
      const NumArray<float> &dx, hid_t my_group, bool with_opendxinfo = true);

  template_noPointer bool Write2DMatrix_withMPI_IO(
      const std::string &ArrayName, const NumMatrix<T, 2> &data, const NumArray<int> &mx_global,
      const NumArray<int> &mx_local, const NumArray<int> &rank_shift, const NumArray<float> &xb,
      const NumArray<float> &dx);
#endif
#endif

  hid_t AddGroup(std::string);
  void CloseGroup(hid_t group);
  hid_t get_defaultGroup();
  /*! Check whether attribute is available in the file*/
  bool doesAttrExist(hid_t, std::string) const;

 private:
  bool AddToEntries();
  bool AddToEntries(hid_t group);
  //	template <typename T>
  //	bool AddAttrToArrSingle(std::string, hid_t &, hid_t, const std::string &, T);
  bool AddDatasetName(const std::string &, hid_t group);
  bool AddDatasetName(const std::string &);
  hid_t hdf5file;
  hid_t group;
  hid_t plist_file_id, plist_dset_id;
  herr_t return_val;
  std::string filename;
  // char* fname;
  int NumEntries;
  int num;
  int rank;
  bool open, use_MPI_IO;
  unsigned MajorNum, MinorNum, ReleaseNum;
};

class Hdf5iStream {
 public:
  Hdf5iStream(std::string, int rank = 0, bool use_MPI_IO = true);
  ~Hdf5iStream();
  bool close();

  template_noPointer bool Read1DMatrix(const std::string &ArrayName, NumMatrix<T, 1> &data);
  template_noPointer bool Read1DMatrix(std::string GroupName, const std::string &ArrayName,
                                       NumMatrix<T, 1> &data);

  template_noPointer bool Read2DMatrix(const std::string &ArrayName, NumMatrix<T, 2> &);
  template_noPointer bool Read2DFrom3DMatrix(const std::string &ArrayName, NumMatrix<T, 2> &data,
                                             int dir, int posPerp);
  template_noPointer bool Read2DFrom3DMatrix(std::string GroupName, const std::string &ArrayName,
                                             NumMatrix<T, 2> &data, int dir, int posPerp);

  template_noPointer bool Read3DMatrix(const std::string &ArrayName, NumMatrix<T, 3> &data,
                                       bool load_size = false);
  template_noPointer int Read3DMatrix(std::string GroupName, const std::string &ArrayName,
                                      NumMatrix<T, 3> &data, bool load_size = false);
  template_noPointer bool Read3DMatrix(const std::string &ArrayName, NumMatrix<T, 3> &data,
                                       double *xb, double *dx);
  template_noPointer bool Read3DMatrix(std::string GroupName, const std::string &ArrayName,
                                       NumMatrix<T, 3> &data, double *xb, double *dx);

#ifdef parallel
  /*!
   * Read data from 3D double file using MPIIO
   * */
  template_noPointer int Read3DMatrix_parallel(std::string GroupName, const std::string &ArrayName,
                                               const NumArray<int> &mx_local,
                                               const NumArray<int> &rank_shift, NumMatrix<T, 3> &);
  /*!
   * Same as above, but for generic group
   * */
  template_noPointer int Read3DMatrix_parallel(const std::string &ArrayName,
                                               const NumArray<int> &mx_local,
                                               const NumArray<int> &rank_shift, NumMatrix<T, 3> &);
#endif

  template_noPointer T ReadPointFrom3DMatrix(const std::string &ArrayName, int ix, int iy, int iz);
  template_noPointer int ReadNumArray(const std::string &ArrayName, NumArray<T> &data);
  template_noPointer int ReadNumArray(const std::string &ArrayName, NumArray<T> &data,
                                      hid_t my_group);
  template_noPointer bool ReadGlobalAttr(const std::string &AttrName, T &AttrData);
  template_noPointer bool ReadGlobalAttr(hid_t my_group, const std::string &AttrName, T &AttrData);
  template_noPointer int ReadGlobalAttr(std::string GroupName, const std::string &AttrName,
                                        T &AttrData);
  template_noPointer T ReadDatasetAttr(std::string DataSetName, const std::string &AttrName);

  /**
   * Read number of grid pointfs from hdf5 dataset
   * */
  void getSize(const std::string &ArrayName, int mx[], int dim);
  void getSize(std::string GroupName, const std::string &ArrayName, int Nx[], int dim);

  bool ReadDatasetGrid(std::string DataSetName, NumArray<float> &xPos, NumArray<float> &yPos,
                       NumArray<float> &zPos);

  int GetDatasetDimension(std::string DataSetName);
  NumArray<int> GetDatasetExtent(std::string DataSetName);
  std::string GetDatasetName(int);
  std::string GetDatasetName(std::string, int);
  // bool doesAttrExist(const char* name ) const;

  bool doesAttrExist(std::string name) const;
  bool doesAttrExist(std::string name, hid_t group) const;
  bool doesGroupExist(std::string name) const;
  bool doesSubGroupExist(std::string name) const;

  bool doesDsetExist(std::string GroupName, std::string DsetName) const;
  bool doesDsetExist(std::string DsetName) const;

  hid_t get_defaultGroup();
  std::string get_defaultGroupName();

#ifdef parallel
#if (HDF_PARALLEL_IO == CRONOS_ON)
  template_noPointer bool Read3DMatrix_withMPI_IO(const std::string &ArrayName,
                                                  NumMatrix<T, 3> &data, NumArray<int> &mx_local,
                                                  NumArray<int> &rank_shift);
  template_noPointer bool Read3DMatrix_withMPI_IO(const std::string &ArrayName,
                                                  NumMatrix<T, 3> &data, NumArray<int> &mx_local,
                                                  NumArray<int> &rank_shift, std::string GroupName);
  template_noPointer bool Read3DMatrix_withMPI_IO(const std::string &ArrayName,
                                                  NumMatrix<T, 3> &data, NumArray<int> &mx_local,
                                                  NumArray<int> &rank_shift, hid_t my_group);
#endif
#endif

 protected:
  hid_t OpenGroup(std::string);
  void CloseGroup(hid_t);

 private:
  hid_t hdf5file;
  hid_t group;
  hid_t plist_file_id, plist_dset_id;
  std::string groupname;
  herr_t return_val;
  unsigned MajorNum, MinorNum, ReleaseNum;
  bool use_MPI_IO;
};

inline bool file_exists(const char *name) {
  std::ifstream f(name);
  return f.good();
}

inline bool file_exists(const std::string name) {
  std::ifstream f(name);
  return f.good();
}

#endif
