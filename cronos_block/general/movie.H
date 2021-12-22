#ifndef CRONOS_MOVIE_H
#define CRONOS_MOVIE_H

#include <stdlib.h>
#include <map>
#include "grid.H"
#include "matrix.H"

/**
 * @brief activates the parallel writing of the mov hdf5
 * @detail parallel io was found to be slower as the rank 0 writing strategy for small datasets
 * I'll deactivate it for now
 */
#define MOV_HDF_PARALLEL_IO CRONOS_OFF

class MovieSlice {
 public:
  MovieSlice(const Grid &grid, const std::string &fieldName, const std::string &fileName);

  void write(const Grid &grid, const NumMatrix<double, 3> &data, const int dirPerp,
             const int globalSliceIndex);

 private:
  /** Computes the direction indices of the parallel to a slicing plane */
  static std::pair<int, int> computeSliceDirs(const int dirPerp);

  /** Returns the slice index in terms of the local grid */
  static int getLocalGridIndex(const Grid &grid, const int dirPerp, const int globalSliceIndex);

  /** Checks if a given slice buffer was already initialized and does so if not */
  void checkAndInitSliceBuffer(const Grid &grid);

  /** Output buffer in order of dirPerp: yz, xz, xy */
  static std::array<NumMatrix<double, 2>, 3> sliceBuffer;

  /** set by setSlice */
  int dirPerp;
  std::pair<int, int> sliceDirs;
  bool overlapsSlice;
  int localSliceIndex;

  /* Set during construction */

  const std::string fieldName;
  std::unique_ptr<Hdf5Stream> h5;

  /* incrementing dataset counter */
  int frameNumber;

  /** sets slice, first to be called in write */
  void setSlice(const Grid &grid, const int dirPerp, const int globalSliceIndex);

  /** Checks if the local rank overlaps with the slicing plane */
  bool checkSliceOverlap(const Grid &grid) const;

  /** Transfers the data slice to the local buffer */
  void fillSliceBuffer(const Grid &grid, const NumMatrix<double, 3> &data);

  /** Returns the name of the dataset to be written */
  std::string getDatasetName(const int outputCount) const;

  /** Writes the hdf5 output of a given slice */
  bool writeHdf5(const Grid &grid, const std::string datasetName,
                 const NumMatrix<double, 2> &sliceData);

#if (defined(parallel) && (MOV_HDF_PARALLEL_IO == CRONOS_OFF))
  /** HDF5 Output buffer for the global slice (rank 0 writing strategy)
   * Only allocated on rank 0
   * subslice data from other ranks is merged
   */
  static std::array<NumMatrix<double, 2>, 3> globalSliceBuffer;

  /** initialize global slice buffer */
  void checkAndInitGlobalSliceBuffer(const Grid &grid);

  /** write a local slice (patch) to the global buffer */
  static void writeLocalToGlobalSlice(const NumMatrix<double, 2> &localSlice,
                                      const std::array<int, 2> &offset,
                                      const std::array<int, 2> &width,
                                      NumMatrix<double, 2> &globalSlice);

  /** Send local slices (patch) to rank 0 */
  void SendDataToMaster(const Grid &grid, NumMatrix<double, 2> &sliceBuffer) const;

  /** recieve local patches */
  void RecvDataFromSlave(const Grid &grid, NumMatrix<double, 2> &globalSliceBuffer) const;
#endif
};

class Movie {
 public:
  Movie();

  /** write standard slices xy, yz, xz  */
  void writeStdSlices(const Grid &grid, const NumMatrix<double, 3> &data,
                      const std::string &fieldName);

  /** write a given slice and initialize the slice if needed  */
  void write(const Grid &grid, const NumMatrix<double, 3> &data, const int dirPerp,
             const int globalSliceIndex, const std::string &fieldName, const std::string &suffix);

 private:
  std::map<std::string, MovieSlice> slices;

  /** Builds the filename */
  static std::string getFilename(const std::string &fieldName, const std::string &suffix);

  /** initialize a given slice */
  void initSlice(const Grid &grid, const std::string &fieldName, const std::string &fileName);

  /** checks if a slice was initialized */
  bool checkInitSlice(const std::string &filename) const;
};

#endif
