#include "movie.H"

#include <sys/stat.h>
#ifdef _MSC_VER
#include "timewrapper.H"
#else
#include <utime.h>
#endif
#include <map>
#include <utility>
#include <filesystem>
#include "CException.H"
#include "Hdf5File_cbase.H"
#include "data.H"

MovieSlice::MovieSlice(const Grid &grid, const std::string &fieldName, const string &fileName)
    : fieldName(fieldName) {
#if (defined(parallel) && (MOV_HDF_PARALLEL_IO == CRONOS_ON))
  h5 =
      std::unique_ptr<Hdf5Stream>(new Hdf5Stream(fileName, 10000, grid.rank, true, MPI_COMM_WORLD));
#else
  if (grid.rank == 0) {
    h5 = std::unique_ptr<Hdf5Stream>(new Hdf5Stream(fileName, 10000));
  }
#endif
  if (h5 != nullptr) h5->close();
  frameNumber = 0;
}

std::array<NumMatrix<double, 2>, 3> MovieSlice::sliceBuffer;

void MovieSlice::setSlice(const Grid &grid, const int dirPerp, const int globalSliceIndex) {
  this->dirPerp = dirPerp;

  sliceDirs = computeSliceDirs(dirPerp);
  localSliceIndex = getLocalGridIndex(grid, dirPerp, globalSliceIndex);
  overlapsSlice = checkSliceOverlap(grid);
}

void MovieSlice::write(const Grid &grid, const NumMatrix<double, 3> &data, const int dirPerp,
                       const int globalSliceIndex) {
  setSlice(grid, dirPerp, globalSliceIndex);

  if (overlapsSlice) {
    checkAndInitSliceBuffer(grid);
    fillSliceBuffer(grid, data);
  }

  const std::string datasetName = getDatasetName(frameNumber);

#if (defined(parallel) && (MOV_HDF_PARALLEL_IO == CRONOS_OFF))
  if (grid.rank == 0) {
    checkAndInitGlobalSliceBuffer(grid);
    RecvDataFromSlave(grid, globalSliceBuffer[dirPerp]);

    h5->Reopen();
    writeHdf5(grid, datasetName, globalSliceBuffer[dirPerp]);
    h5->close();
  } else {
    SendDataToMaster(grid, sliceBuffer[dirPerp]);
  }
#else
  h5->Reopen();
  writeHdf5(grid, datasetName, sliceBuffer[dirPerp]);
  h5->close();
#endif

  ++frameNumber;
}

std::pair<int, int> MovieSlice::computeSliceDirs(const int dirPerp) {
  assert(dirPerp >= 0);
  assert(dirPerp <= 2);
  switch (dirPerp) {
    case 0:
      return std::make_pair(1, 2);
    case 1:
      return std::make_pair(0, 2);
    case 2:
      return std::make_pair(0, 1);
    default:
      throw CException("Wrong direction! " + to_string(dirPerp));
  }
}

int MovieSlice::getLocalGridIndex(const Grid &grid, const int dirPerp, const int globalSliceIndex) {
  return globalSliceIndex - grid.get_RankShift(dirPerp);
}

bool MovieSlice::checkSliceOverlap(const Grid &grid) const {
  return (0 <= localSliceIndex) && (localSliceIndex <= grid.mx[dirPerp]);
}

void MovieSlice::checkAndInitSliceBuffer(const Grid &grid) {
  if (sliceBuffer[dirPerp].getSize() == 0) {
    sliceBuffer[dirPerp].resize(Index::set(0, 0),
                                Index::set(grid.mx[sliceDirs.first], grid.mx[sliceDirs.second]));
  }
}

void MovieSlice::fillSliceBuffer(const Grid &grid, const NumMatrix<double, 3> &data) {
  assert(localSliceIndex >= 0);
  assert(localSliceIndex <= grid.mx[dirPerp]);
  assert(sliceBuffer[dirPerp].getSize() ==
         (grid.mx[sliceDirs.first] + 1) * (grid.mx[sliceDirs.second] + 1));

  int idx[3];
  idx[dirPerp] = localSliceIndex;
  int &iDir1 = idx[sliceDirs.first];
  int &iDir2 = idx[sliceDirs.second];

  for (iDir2 = 0; iDir2 <= grid.mx[sliceDirs.second]; ++iDir2) {
    for (iDir1 = 0; iDir1 <= grid.mx[sliceDirs.first]; ++iDir1) {
      sliceBuffer[dirPerp](iDir1, iDir2) = data[idx];
    }
  }
}

std::string MovieSlice::getDatasetName(const int outputCount) const {
  return fieldName + ToString(outputCount, 5);
}

bool MovieSlice::writeHdf5(const Grid &grid, const std::string datasetName,
                           const NumMatrix<double, 2> &sliceData) {
  bool ret = true;

  auto xbSliceGlobal =
      NumArray<double>::set(grid.global_xb[sliceDirs.first], grid.global_xb[sliceDirs.second]);
  auto dxSliceGlobal = NumArray<double>::set(grid.dx[sliceDirs.first], grid.dx[sliceDirs.second]);

#if (defined(parallel) && (MOV_HDF_PARALLEL_IO == CRONOS_ON))
  NumArray<int> mxSliceGlobal =
      Index::set(grid.global_mx[sliceDirs.first], grid.global_mx[sliceDirs.second]);

  NumArray<int> mxSliceLocal;
  if (overlapsSlice) {
    mxSliceLocal = Index::set(grid.mx[sliceDirs.first], grid.mx[sliceDirs.second]);
  } else {
    mxSliceLocal = Index::set(-1, -1);  // To indicate that this rank writes no data
  }

  NumArray<int> rankShift =
      Index::set(grid.get_RankShift(sliceDirs.first), grid.get_RankShift(sliceDirs.second));

  ret =
      h5->Write2DMatrix_withMPI_IO(datasetName, sliceData, mxSliceGlobal, mxSliceLocal, rankShift,
                                   xbSliceGlobal.convert<float>(), dxSliceGlobal.convert<float>());
#else

  ret = h5->Write2DMatrix(datasetName, sliceData, xbSliceGlobal.data(), dxSliceGlobal.data(), true);
#endif
  h5->AddAttributeToArray(datasetName, "num", frameNumber);

  // Write also time if possible
  if (const Data *gdata = dynamic_cast<const Data *>(&grid)) {
    h5->AddAttributeToArray(datasetName, "time", gdata->time);
  }

  return ret;
}

#if (defined(parallel) && (MOV_HDF_PARALLEL_IO == CRONOS_OFF))
std::array<NumMatrix<double, 2>, 3> MovieSlice::globalSliceBuffer;

void MovieSlice::checkAndInitGlobalSliceBuffer(const Grid &grid) {
  assert(grid.rank == 0);

  if (globalSliceBuffer[dirPerp].getSize() == 0) {
    globalSliceBuffer[dirPerp].resize(
        Index::set(0, 0),
        Index::set(grid.global_mx[sliceDirs.first], grid.global_mx[sliceDirs.second]));
  }
}

void MovieSlice::writeLocalToGlobalSlice(const NumMatrix<double, 2> &localSlice,
                                         const std::array<int, 2> &offset,
                                         const std::array<int, 2> &width,
                                         NumMatrix<double, 2> &globalSlice) {
  for (int iDir2 = 0; iDir2 < width[1]; ++iDir2) {
    for (int iDir1 = 0; iDir1 < width[0]; ++iDir1) {
      globalSlice(iDir1 + offset[0], iDir2 + offset[1]) = localSlice(iDir1, iDir2);
    }
  }
}

enum CommTag { SendFlag, GridStructure, Data };

void MovieSlice::SendDataToMaster(const Grid &grid, NumMatrix<double, 2> &sliceBuffer) const {
  assert(grid.rank != 0);

  if (overlapsSlice) {
    int sndFlagBuffer = 1;
    MPI_Send(&sndFlagBuffer, 1, MPI_INT, 0, CommTag::SendFlag, MPI_COMM_WORLD);

    std::array<int, 4> gridStructure = {
        grid.get_RankShift(sliceDirs.first), grid.get_RankShift(sliceDirs.second),
        grid.get_RankWidth(sliceDirs.first), grid.get_RankWidth(sliceDirs.second)};
    MPI_Send(gridStructure.data(), 4, MPI_INT, 0, CommTag::GridStructure, MPI_COMM_WORLD);

    MPI_Send((double *)sliceBuffer, sliceBuffer.getSize(), MPI_DOUBLE, 0, CommTag::Data,
             MPI_COMM_WORLD);
  } else {
    int sndFlagBuffer = 0;
    MPI_Send(&sndFlagBuffer, 1, MPI_INT, 0, CommTag::SendFlag, MPI_COMM_WORLD);
  }
}

void MovieSlice::RecvDataFromSlave(const Grid &grid,
                                   NumMatrix<double, 2> &globalSliceBuffer) const {
  assert(grid.rank == 0);
  assert(globalSliceBuffer.getSize() ==
         (grid.global_mx[sliceDirs.first] + 1) * (grid.global_mx[sliceDirs.second] + 1));

  if (overlapsSlice) {
    std::array<int, 2> offset = {grid.get_RankShift(sliceDirs.first),
                                 grid.get_RankShift(sliceDirs.second)};
    std::array<int, 2> width = {grid.get_RankWidth(sliceDirs.first),
                                grid.get_RankWidth(sliceDirs.second)};
    writeLocalToGlobalSlice(sliceBuffer[dirPerp], offset, width, globalSliceBuffer);
  }

  NumMatrix<double, 2> rcvBuffer;
  for (int sendingRank = 1; sendingRank < grid.ntasks; sendingRank++) {
    int rcvFlag;

    MPI_Recv(&rcvFlag, 1, MPI_INT, sendingRank, CommTag::SendFlag, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);

    if (rcvFlag) {
      std::array<int, 4> gridStructure;
      MPI_Recv(gridStructure.data(), 4, MPI_INT, sendingRank, CommTag::GridStructure,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      if (rcvBuffer.getDims(0) != gridStructure[2] || rcvBuffer.getDims(1) != gridStructure[3]) {
        rcvBuffer.resize(&gridStructure[2]);
      }

      MPI_Recv((double *)rcvBuffer, rcvBuffer.getSize(), MPI_DOUBLE, sendingRank, CommTag::Data,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      writeLocalToGlobalSlice(rcvBuffer, {gridStructure[0], gridStructure[1]},
                              {gridStructure[2], gridStructure[3]}, globalSliceBuffer);
    }
  }
}

#endif

Movie::Movie(){};

std::string Movie::getFilename(const std::string &fieldName, const std::string &suffix) {
  char dirname[255];
  std::ostringstream dummy;

  sprintf(dirname, "%s/%s_mov", getenv("poub"), getenv("pname"));
  //mkdir(dirname, 511);
  filesystem::create_directory(dirname);
  filesystem::permissions(dirname, filesystem::perms::owner_read | filesystem::perms::owner_write | filesystem::perms::owner_exec | filesystem::perms::group_exec | filesystem::perms::others_exec);
  //utime(dirname, NULL);
  auto ftime = filesystem::last_write_time(dirname);
  filesystem::file_time_type currentTime(decltype(ftime)::clock::now());
  filesystem::last_write_time(dirname, currentTime);

  dummy << dirname << "/" << getenv("pname") << "_" << fieldName << "_" << suffix << ".h5";
  return dummy.str();
}

void Movie::writeStdSlices(const Grid &grid, const NumMatrix<double, 3> &data,
                           const std::string &fieldName) {
  write(grid, data, 0, grid.global_mx[0] / 2, fieldName, "movyz");
  write(grid, data, 1, grid.global_mx[1] / 2, fieldName, "movxz");
  write(grid, data, 2, grid.global_mx[2] / 2, fieldName, "movxy");
}

void Movie::write(const Grid &grid, const NumMatrix<double, 3> &data, const int dirPerp,
                  const int globalSliceIndex, const std::string &fieldName,
                  const std::string &suffix) {
  const std::string fileName = getFilename(fieldName, suffix);
  if (!checkInitSlice(fileName)) {
    initSlice(grid, fieldName, fileName);
  }

  slices.at(fileName).write(grid, data, dirPerp, globalSliceIndex);
}

void Movie::initSlice(const Grid &grid, const std::string &fieldName, const std::string &fileName) {
  slices.insert(std::make_pair(fileName, MovieSlice(grid, fieldName, fileName)));
}

bool Movie::checkInitSlice(const std::string &fileName) const {
  return !(slices.find(fileName) == slices.end());
}
