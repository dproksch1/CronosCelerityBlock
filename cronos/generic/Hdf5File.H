#ifndef CRONOS_HDF5FILE_H
#define CRONOS_HDF5FILE_H = 1



#include "H5Cpp.h"
#include "matrix.H"
class Hdf5Stream{
public:
	Hdf5Stream(std::string, int, int rank=0);
	~Hdf5Stream();
	bool close();
	void Reopen();
	bool AddGlobalAttr(std::string AttrName, double);
	bool ChangeGlobalAttr(std::string AttrName, double);
	bool AddGlobalAttr(std::string AttrName, float);
	bool ChangeGlobalAttr(std::string AttrName, float);
	bool AddGlobalAttr(std::string AttrName, int);
	bool ChangeGlobalAttr(std::string AttrName, int);
	bool AddGlobalAttr(std::string AttrName, std::string AttrData);
	// bool AddGlobalAttr(std::string AttrName, std::string *AttrData);
	bool AddGlobalAttr(std::string AttrName, double *AttrData, int);
	bool AddGlobalAttr(std::string AttrName, float *AttrData, int);
	bool AddGlobalAttr(std::string AttrName, unsigned long *AttrData, int);
	bool AddGlobalAttr(std::string AttrName, int *AttrData, int);
	bool AddGlobalAttr(std::string AttrName, std::string *AttrData, int num);
	bool AddAttributeToArray(std::string, const std::string &, int);
	bool AddAttributeToArray(std::string, const std::string &, float);
	bool AddAttributeToArray(std::string, const std::string &, double);

	bool Write1DMatrix(std::string ArrayName, NumMatrix<double,1> &data,
	                   double xb, double dx, int);
	bool Write1DMatrix(std::string ArrayName, NumMatrix<float,1> &data);
	bool Write3DMatrix(std::string ArrayName, NumMatrix<double,3> &data,
	                   double *xb, double *dx);
	bool Write3DVecMatrix(std::string ArrayName, NumMatrix<float,3> &data_x,
	                      NumMatrix<float,3> &data_y, NumMatrix<float,3> &data_z);
	std::string Write2DMatrix(std::string ArrayName, NumMatrix<double,2> &data,
	                          double *xb, double *dx, int);
	bool Write3DMatrix(std::string ArrayName, NumMatrix<float,3> &data,
	                   double *xb, double *dx, bool with_opendxinfo=true);
	bool Write3DMatrix(std::string ArrayName, NumMatrix<float,3> &data);
	bool Write3DMatrixSwap(std::string ArrayName, NumMatrix<double,3> &data,
	                       double *xb, double *dx);
	bool Write2DMatrix(std::string ArrayName, NumMatrix<float,2> &data,
	                   double *xb, double *dx, bool with_opendxinfo=true);
	bool WriteArray(int *data, int max);
	bool WriteArray(std::string ,int *data, int max);
	bool WriteArray(float *data, int max);
	bool WriteArray(std::string ,float *data, int max);
	bool WriteArray(double *data, int max);
	bool WriteArray(std::string ,double *data, int max);
	bool WriteNDArray(std::string ,float *data, int mx[], int dim);
private:
	template <typename T>
	bool AddAttrToArrSingle(std::string, H5::AtomType &, const std::string &, T);
	bool AddDatasetName(std::string &);
	H5::H5File* hdf5file;
	H5::Group* group;
	std::string filename;
	int NumEntries;
	int num;
	bool open;
	unsigned MajorNum, MinorNum, ReleaseNum;
};



class Hdf5iStream{
public:
	Hdf5iStream(std::string, int rank=0);
	~Hdf5iStream();
	bool close();
	bool Read3DMatrix(std::string, NumMatrix<float,3> &);
	bool Read3DMatrix(std::string ArrayName, NumMatrix<float,3> &data,
	                  double *xb, double *dx);
	bool Read3DMatrix(std::string, NumMatrix<double,3> &);
	void getSize(std::string, int [], int);
	bool ReadGlobalAttr(std::string , double &);
	bool ReadGlobalAttr(std::string , float &);
	bool ReadGlobalAttr(std::string , unsigned long &);
	bool ReadGlobalAttr(std::string , int &);
	std::string GetDatasetName(int);
	bool doesAttrExist(const char* name ) const;
private:
	H5::H5File* hdf5file;
	H5::Group* group;
	int NumEntries;
	unsigned MajorNum, MinorNum, ReleaseNum;
};

#endif
