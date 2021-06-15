/*
  Deprecated version of hdf5Stream class
  --> functionality now included in Hdf5File_cbase.C
  --> Current file is to be deleted from repo in the future

 */

#include <iostream>
#include <stdlib.h>
#include "Hdf5File.H"

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif
using namespace std;

#define DELETE(x) { if (x != NULL) delete x; x = NULL; }

Hdf5Stream::Hdf5Stream(string filename, int NumEntries, int rank)
{
	this->filename = filename;
#ifdef AIX
	int len = filename.size();
	char* fname = new char[len+1];
	filename.copy(fname,len);
	fname[len] = '\0';
	hdf5file = new H5File(fname, H5F_ACC_TRUNC );
	delete [] fname;
#else
	hdf5file = new H5File(filename, H5F_ACC_TRUNC );
#endif
	// All Data to be written into group data (this allows for global
	// attributes)
	group    = new Group(hdf5file->createGroup( "/Data" ));

	IntType datatype( PredType::NATIVE_INT );
	hsize_t DimsInfo = 1;
	DataSpace infospace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute("Entries", datatype, infospace));
	Info.write( datatype, &NumEntries);
  
	this->NumEntries = NumEntries;
	this->num = 0;
	this->open = true;

	// Getting the version 
	H5::H5Library::getLibVersion(MajorNum, MinorNum, ReleaseNum);

}



void Hdf5Stream::Reopen()
{
	if(!open) {
	  DELETE(hdf5file);
	  DELETE(group);
		hdf5file = new H5File(filename, H5F_ACC_RDWR );
		group = new Group( hdf5file->openGroup( "Data" ));
	}

	this->open = true;
}



Hdf5iStream::Hdf5iStream(string filename, int rank)
{
	if(rank == 0){
		cout << " Opening file: " << filename;
	}
	hdf5file = new H5File(filename, H5F_ACC_RDONLY );
	group = new Group( hdf5file->openGroup( "Data" ));
	if(rank == 0){
		cout << " ...done " << endl;
	}

	int numh5;
	Attribute Info = group->openAttribute( "Entries");
	Info.read(PredType::NATIVE_INT, &numh5);

	this->NumEntries = numh5;
	// Getting the version 
	H5::H5Library::getLibVersion(MajorNum, MinorNum, ReleaseNum);

}




bool Hdf5Stream::AddGlobalAttr(string AttrName, double AttrData)
{
	FloatType datatype( PredType::NATIVE_DOUBLE );
	hsize_t DimsInfo = 1;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, &AttrData);
	return true;
}


bool Hdf5Stream::ChangeGlobalAttr(string AttrName, double AttrData)
{
	Attribute Info = group->openAttribute(AttrName);
	Info.write( PredType::NATIVE_DOUBLE, &AttrData);
	return true;
}


bool Hdf5iStream::ReadGlobalAttr(string AttrName, double &AttrData)
{
	Attribute Info = group->openAttribute(AttrName);
	Info.read(PredType::NATIVE_DOUBLE, &AttrData);
	return true;
}


bool Hdf5Stream::AddGlobalAttr(string AttrName, float AttrData)
{
	FloatType datatype( PredType::NATIVE_FLOAT );
	hsize_t DimsInfo = 1;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, &AttrData);
	return true;
}

bool Hdf5Stream::ChangeGlobalAttr(string AttrName, float AttrData)
{
	Attribute Info = group->openAttribute(AttrName);
	Info.write( PredType::NATIVE_FLOAT, &AttrData);
	return true;
}

bool Hdf5iStream::ReadGlobalAttr(string AttrName, float &AttrData)
{
	Attribute Info = group->openAttribute(AttrName);
	Info.read(PredType::NATIVE_FLOAT, &AttrData);
	return true;
}


bool Hdf5Stream::AddGlobalAttr(string AttrName, int AttrData)
{
	IntType datatype( PredType::NATIVE_INT );
	hsize_t DimsInfo = 1;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, &AttrData);
	return true;
}


bool Hdf5Stream::ChangeGlobalAttr(string AttrName, int AttrData)
{
	Attribute Info = group->openAttribute(AttrName);
	Info.write( PredType::NATIVE_INT, &AttrData);
	return true;
}

bool Hdf5iStream::ReadGlobalAttr(string AttrName, int &AttrData)
{

	Attribute Info = group->openAttribute(AttrName);
	
	Info.read(PredType::NATIVE_INT, &AttrData);
	return true;
}

bool Hdf5iStream::ReadGlobalAttr(string AttrName, unsigned long &AttrData)
{
	Attribute Info = group->openAttribute(AttrName);
	Info.read(PredType::NATIVE_ULONG, &AttrData);
	return true;
}


bool Hdf5Stream::AddGlobalAttr(string AttrName, string AttrData)
{
	// Set datatype to string
	StrType datatype( PredType::C_S1);
	if(ReleaseNum > 5 || MinorNum >= 8) {
		datatype.setSize( H5T_VARIABLE );
	} else if (ReleaseNum == 5) {
		// Get length of string
		hsize_t StrLen(0);
		for(int i=0; i<num; ++i) {
			StrLen = AttrData.size();
		}
		datatype.setSize(StrLen+1);
	}

	hsize_t DimsInfo = 1;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, &AttrData);
	return true;
}


bool Hdf5Stream::AddGlobalAttr(string AttrName, double *AttrData, int num)
{
	FloatType datatype( PredType::NATIVE_DOUBLE );
	hsize_t DimsInfo = num;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, AttrData);
	return true;
}



bool Hdf5Stream::AddGlobalAttr(string AttrName, float *AttrData, int num)
{
	FloatType datatype( PredType::NATIVE_FLOAT );
	hsize_t DimsInfo = num;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, AttrData);
	return true;
}


bool Hdf5Stream::AddGlobalAttr(string AttrName,unsigned long *AttrData, int num)
{
	IntType datatype( PredType::NATIVE_ULONG );
	hsize_t DimsInfo = num;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, AttrData);
	return true;
}


bool Hdf5Stream::AddGlobalAttr(string AttrName, int *AttrData, int num)
{
	IntType datatype( PredType::NATIVE_INT );
	hsize_t DimsInfo = num;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, AttrData);
	return true;
}


bool Hdf5Stream::AddGlobalAttr(string AttrName, string *AttrData, int num)
{
	// Set datatype to string
	StrType datatype( PredType::C_S1);
	// Set length of string according to hdf5 version
	cout << AttrData[0] << endl << AttrData[1] << endl;
	if(ReleaseNum > 5 || MinorNum >= 8) {
		datatype.setSize( H5T_VARIABLE );
	} else if (ReleaseNum == 5) {
		// Get length of string
		hsize_t StrLen(0);
		for(int i=0; i<num; ++i) {
			StrLen = AttrData[i].size();
		}
		datatype.setSize(StrLen+1);
	}

	hsize_t DimsInfo = num;
	DataSpace AttrSpace( 1, &DimsInfo);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, AttrData);
	return true;
}

bool Hdf5Stream::AddDatasetName(string &DatasetName) {

	char cnum[255];
	sprintf(cnum,"%2.2d",this->num);
	string AttrName = "Name_om";
	AttrName += cnum;
 
	hsize_t StrLen = DatasetName.size();
	StrType datatype( PredType::C_S1);
	//  StrType datatype( H5T_C_S1 );

	if(ReleaseNum > 5 || MinorNum >= 8) {
		datatype.setSize( H5T_VARIABLE );
	} else if (ReleaseNum == 5) {
		datatype.setSize(StrLen+1);
	}
	//   datatype.setSize(StrLen);
	//   hsize_t DimsInfo = 1;
	// /
	//   DataSpace AttrSpace( 1, &DimsInfo);
	DataSpace AttrSpace( H5S_SCALAR );
	//   StrType datatype( PredType::C_S1);
	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute Info(group->createAttribute(AttrName, datatype, AttrSpace));
	Info.write( datatype, DatasetName);
	//   string fu;
	//   Info.read(datatype, fu);
	//   StrType datatypeVar( PredType::C_S1 , Info.getStorageSize());
	return true;
}

string Hdf5iStream::GetDatasetName(int num) {
	char cnum[255];
	sprintf(cnum,"%2.2d",num);
	string AttrName = "Name_om";
	AttrName += cnum;

	Attribute Info = group->openAttribute(AttrName);
	string DatasetName;
	StrType datatype( PredType::C_S1 );
	//  if(ReleaseNum == 6) {
	if(ReleaseNum > 5 || MinorNum >= 8) {
		datatype.setSize( H5T_VARIABLE );
	} else if (ReleaseNum == 5) {
		datatype.setSize(Info.getStorageSize());
	}
	Info.read(datatype, DatasetName);
	return DatasetName;
}





bool Hdf5Stream::Write1DMatrix(string ArrayName, NumMatrix<double,1> &data,
                               double Origin, double Delta, int numin)
{
	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.

	   Remarks:
     
	   On reading the hdf data the swapped dimensions have to be taken
	   into account
    
	*/
	int mx = data.getHigh(0) - data.getLow(0) + 1;
  
	num+=1;

	int DIM = 1;
	FloatType datatype( PredType::NATIVE_DOUBLE );
	datatype.setOrder( H5T_ORDER_LE );
  
	hsize_t DimsData = mx;
	DataSpace dataspace( DIM, &DimsData);
 

	// Supplying additional attributes for opendx input

	FloatType datatypefloat( PredType::NATIVE_DOUBLE );
	hsize_t DimsAttr = 1;
	DataSpace attrspace( 1, &DimsAttr );


	DataSet dataset = group->createDataSet( ArrayName, datatype, dataspace );
	// Writing attributes
	Attribute origin = dataset.createAttribute("origin", datatypefloat, attrspace);
	origin.write( PredType::NATIVE_DOUBLE, &Origin);
	Attribute delta = dataset.createAttribute("delta", datatypefloat, attrspace);
	delta.write( datatypefloat, &Delta);
	dataset.write( data, datatype );
	return true;
}



bool Hdf5Stream::Write1DMatrix(string ArrayName, NumMatrix<float,1> &data)
{
	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.

	   Remarks:
     
	*/
	int mx = data.getHigh(0) - data.getLow(0) + 1;
  
	num+=1;

	int DIM = 1;
	FloatType datatype( PredType::NATIVE_FLOAT );
	datatype.setOrder( H5T_ORDER_LE );
  
	// Create dataspace
	hsize_t DimsData = mx;
	DataSpace dataspace( DIM, &DimsData);
 
	// Create dataset
	DataSet dataset = group->createDataSet( ArrayName, datatype, dataspace );

	// Write data
	dataset.write( data, datatype );
	return true;
}



// bool Hdf5Stream::AddAttributeToArray(string ArrayName,
// 				     const string &AttributeName,
// 				     double AttributeData)
// {
//   /*******************************************************
//    *  Routine to add Attribute to Any written Array Data *
//    ******************************************************/

//   ArrayName = "/Data/"+ArrayName;

//   // Reopen written dataset:
//   DataSet dataset = hdf5file->openDataSet( ArrayName );


//   // Preparing Attribute:
//   FloatType datatype( PredType::NATIVE_DOUBLE );
//   hsize_t DimsInfo = 1;
//   DataSpace AttributeSpace( 1, &DimsInfo);
//   datatype.setOrder( H5T_ORDER_LE ); // Little endian
//   Attribute attribute = group->createAttribute( AttributeName, datatype, AttributeSpace);

//   return true;

// }

bool Hdf5Stream::AddAttributeToArray(string ArrayName,
                                     const string &AttributeName,
                                     double AttributeData)
{

	/*******************************************************
	 *  Routine to write double Attribute to Array Data   *
	 ******************************************************/

	// Preparing Attribute:
	FloatType datatype( PredType::NATIVE_DOUBLE );

	return AddAttrToArrSingle(ArrayName, datatype, AttributeName, AttributeData);
  
}


bool Hdf5Stream::AddAttributeToArray(string ArrayName,
                                     const string &AttributeName,
                                     float AttributeData)
{

	/*******************************************************
	 *  Routine to write float Attribute to Array Data   *
	 ******************************************************/

	// Preparing Attribute:
	FloatType datatype( PredType::NATIVE_FLOAT );

	return AddAttrToArrSingle(ArrayName, datatype, AttributeName, AttributeData);
  
}



bool Hdf5Stream::AddAttributeToArray(string ArrayName,
                                     const string &AttributeName,
                                     int AttributeData)
{

	/*******************************************************
	 *  Routine to write integer Attribute to Array Data   *
	 ******************************************************/
	// Preparing Attribute:
	IntType datatype( PredType::NATIVE_INT );

	return AddAttrToArrSingle(ArrayName, datatype, AttributeName, AttributeData);
  
}

template <typename T>
bool Hdf5Stream::AddAttrToArrSingle(string ArrayName,
                                    AtomType &datatype,
                                    const string &AttributeName,
                                    T AttributeData)
{
	/*******************************************************
	 *  Routine to add Attribute to Any written Array Data *
	 *  the corresponding array is identified by its name  *
	 ******************************************************/

	// Reopen written dataset:
	DataSet dataset = group->openDataSet( ArrayName );

	hsize_t DimsInfo = 1;
	DataSpace AttributeSpace( 1, &DimsInfo);

	datatype.setOrder( H5T_ORDER_LE ); // Little endian
	Attribute attribute = dataset.createAttribute( AttributeName, datatype, AttributeSpace);
	attribute.write( datatype, &AttributeData );
	return true;
}



bool Hdf5Stream::Write3DMatrix(string ArrayName, NumMatrix<double,3> &data,
                               double *xb, double *dx)
{
	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.

	   Remarks:
     
	   On reading the hdf data the swapped dimensions have to be taken
	   into account
    
	*/
	int mx[3];
	mx[0]=data.getHigh(2) - data.getLow(2) + 1;
	mx[1]=data.getHigh(1) - data.getLow(1) + 1;
	mx[2]=data.getHigh(0) - data.getLow(0) + 1;
  
	AddDatasetName(ArrayName);

	num+=1;

	const int DIM = 3;
	FloatType datatype( PredType::NATIVE_DOUBLE );
	datatype.setOrder( H5T_ORDER_LE );
  
	hsize_t DimsData[DIM];
	for(int q=0; q<DIM; ++q){
		DimsData[q]  = mx[q];
	}
	DataSpace dataspace( DIM, DimsData);
 

	// Supplying additional attributes for opendx input

	FloatType datatypefloat( PredType::NATIVE_DOUBLE );
	hsize_t DimsAttr = 3;
	DataSpace attrspace( 1, &DimsAttr );

	double Origin[3];
	double Delta[3];
	for(int q=0; q<3; ++q){
		Origin[q] = xb[q];
		Delta[q]  = dx[q];
	}

	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));
	// Writing attributes
	Attribute origin(dataset.createAttribute("origin", datatypefloat, attrspace));
	Attribute delta(dataset.createAttribute("delta", datatypefloat, attrspace));
	origin.write( datatypefloat, Origin);
	delta.write( datatypefloat, Delta);

	dataset.write( data, datatype);

	return true;
}




string Hdf5Stream::Write2DMatrix(string ArrayName, NumMatrix<double,2> &data,
                                 double *xb, double *dx, int numin)
{
	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.

	   Remarks:
     
	   On reading the hdf data the swapped dimensions have to be taken
	   into account
    
	*/
	int mx[2];
	mx[0]=data.getHigh(1) - data.getLow(1) + 1;
	mx[1]=data.getHigh(0) - data.getLow(0) + 1;

	string SaveName = ArrayName;

	char numchar[255];
	sprintf(numchar,"%5.5i",numin);
	ArrayName = ArrayName+numchar;
	SaveName = SaveName+numchar;
	num+=1;
	if(num > NumEntries){
		return false;
	}

	const int DIM = 2;
	FloatType datatype( PredType::NATIVE_DOUBLE );
	datatype.setOrder( H5T_ORDER_LE );
  
	hsize_t DimsData[DIM];
	for(int q=0; q<DIM; ++q){
		DimsData[q]  = mx[q];
	}
	DataSpace dataspace( DIM, DimsData);
 

	// Supplying additional attributes for opendx input

	FloatType datatypefloat( PredType::NATIVE_DOUBLE );
	hsize_t DimsAttr = 2;
	DataSpace attrspace( 1, &DimsAttr );

	double Origin[2];
	double Delta[2];
	for(int q=0; q<DIM; ++q){
		Origin[q] = xb[q];
		Delta[q]  = dx[q];
	}


	// Supplying frame number for movie in opendx:
	IntType datatypenum( PredType::NATIVE_INT );
	hsize_t DimsInfo = 1;
	DataSpace numspace( 1, &DimsInfo);
	datatypenum.setOrder( H5T_ORDER_LE ); // Little endian

	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));
	// Writing attributes
	Attribute origin(dataset.createAttribute("origin", datatypefloat, attrspace));
	Attribute delta(dataset.createAttribute("delta", datatypefloat, attrspace));

	origin.write( datatypefloat, Origin);
	delta.write( datatypefloat, Delta);

	Attribute number(dataset.createAttribute("num", datatypenum, numspace));
	number.write(datatypenum, &num);

	dataset.write( data, datatype );

	return SaveName;
}



bool Hdf5Stream::Write3DVecMatrix(string ArrayName,
                                  NumMatrix<float,3> &data_x,
                                  NumMatrix<float,3> &data_y,
                                  NumMatrix<float,3> &data_z)
{
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
	mx[0]=data_x.getHigh(2) - data_x.getLow(2) + 1;
	mx[1]=data_x.getHigh(1) - data_x.getLow(1) + 1;
	mx[2]=data_x.getHigh(0) - data_x.getLow(0) + 1;
  
	AddDatasetName(ArrayName);

	num+=1;

	const int DIM = 3;

	/* 
     * Define array datatype for the data in the file.
     - 1 component 
     - 3 entries
     */ 
	hsize_t ArrayExtent[1];
	ArrayExtent[0] = 3;
    ArrayType datatype( PredType::NATIVE_FLOAT, 1, ArrayExtent);
 

	// FloatType datatype( PredType::NATIVE_DOUBLE );
	// datatype.setOrder( H5T_ORDER_LE );
  
	hsize_t DimsData[DIM];
	for(int q=0; q<DIM; ++q){
		DimsData[q]  = mx[q];
	}
	DataSpace dataspace( DIM, DimsData);
 
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
				data[i][j][k][0] = data_x(k+ib[2],j+ib[1],i+ib[0]);
				data[i][j][k][1] = data_y(k+ib[2],j+ib[1],i+ib[0]);
				data[i][j][k][2] = data_z(k+ib[2],j+ib[1],i+ib[0]);
			}
		}
	}

	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));
	dataset.write( data.data(), datatype);
	return true;
}






bool Hdf5iStream::Read3DMatrix(string ArrayName, NumMatrix<double,3> &data)
{
	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.

	   Remarks:
     
	   On reading the hdf data the swapped dimensions have to be taken
	   into account
    
	*/

	const int DIM = 3;
	int mx[3];
	mx[0] = data.getHigh(2)-data.getLow(2)+1;
	mx[1] = data.getHigh(1)-data.getLow(1)+1;
	mx[2] = data.getHigh(0)-data.getLow(0)+1;

	DataSet dataset = group->openDataSet( ArrayName );
	DataSpace dataspace = dataset.getSpace();

	int dimhdf = dataspace.getSimpleExtentNdims();
	if(dimhdf != DIM){
		cerr << " Wrong dimensionality of input data: " << endl;
		cerr << dimhdf << " " << DIM << endl;
		exit(-2);
	}
	hsize_t dims_out[DIM];
	int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
	if(ndims != DIM) {
		cerr << " Wrong number of dimensions " << ndims << " - " << DIM << endl;
		exit(-22);
	}
	for(int i=0; i<DIM; ++i){
		if((int(dims_out[i])) != mx[i]){
			cerr << " Wrong size of dimension " << i << ":" << endl;
			cerr << int(dims_out[i]) << " " << mx[i] << endl;
		}
	}

	dataset.read(data, PredType::NATIVE_DOUBLE );
	return true;
}


void Hdf5iStream::getSize(string ArrayName, int mx[], int DIM)
{
  
	DataSet dataset = group->openDataSet( ArrayName );
	DataSpace dataspace = dataset.getSpace();

	int dimhdf = dataspace.getSimpleExtentNdims();
	if(dimhdf != DIM){
		cerr << " Wrong dimensionality of input data: " << endl;
		cerr << dimhdf << " " << DIM << endl;
		exit(-2);
	}

	std::vector<hsize_t> dims_out(DIM);
	dataspace.getSimpleExtentDims( dims_out.data(), NULL);
	for(int i=0; i<DIM; ++i){
		mx[DIM-i-1] = dims_out[i];
	}
}



bool Hdf5Stream::Write3DMatrix(string ArrayName, NumMatrix<float,3> &data) {
	double dummy[3];
	return Write3DMatrix(ArrayName, data, dummy, dummy, false);
}


bool Hdf5Stream::Write3DMatrix(string ArrayName, NumMatrix<float,3> &data,
                               double *xb, double *dx, bool with_opendxinfo)
{
	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.

	   Remarks:
     
	   On reading the hdf data the swapped dimensions have to be taken
	   into account
    
	*/
	int mx[3];
	mx[0]=data.getHigh(2) - data.getLow(2) + 1;
	mx[1]=data.getHigh(1) - data.getLow(1) + 1;
	mx[2]=data.getHigh(0) - data.getLow(0) + 1;
  
	AddDatasetName(ArrayName);

	num+=1;

	const int DIM = 3;
	FloatType datatype( PredType::NATIVE_FLOAT );
	datatype.setOrder( H5T_ORDER_LE );
  
	hsize_t DimsData[DIM];
	for(int q=0; q<DIM; ++q){
		DimsData[q]  = mx[q];
	}
	DataSpace dataspace( DIM, DimsData);
 
	// Supplying additional attributes for opendx input
	FloatType datatypefloat( PredType::NATIVE_FLOAT );
	hsize_t DimsAttr = 3;
	DataSpace attrspace( 1, &DimsAttr );
	float Origin[3];
	float Delta[3];
	if( with_opendxinfo ) {
		for(int q=0; q<3; ++q){
			Origin[q] = float(xb[q]);
			Delta[q]  = float(dx[q]);
		}
	}

	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));

	// Writing attributes (if necessary)
	if( with_opendxinfo ) {
		Attribute origin(dataset.createAttribute("origin", datatypefloat, attrspace));
		Attribute delta(dataset.createAttribute("delta", datatypefloat, attrspace));
		origin.write( datatypefloat, Origin);
		delta.write( datatypefloat, Delta);
	}
	dataset.write( data, datatype );

	return true;
}



bool Hdf5iStream::Read3DMatrix(string ArrayName, NumMatrix<float,3> &data)
{
	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.

	   Remarks:
     
	   On reading the hdf data the swapped dimensions have to be taken
	   into account
    
	*/

	const int DIM = 3;
	int lbound[3], ubound[3];

	DataSet dataset = group->openDataSet( ArrayName );
	DataSpace dataspace = dataset.getSpace();

	int dimhdf = dataspace.getSimpleExtentNdims();
	if(dimhdf != DIM){
		cerr << " Wrong dimensionality of input data: " << endl;
		cerr << dimhdf << " " << DIM << endl;
		exit(-2);
	}
	hsize_t dims_out[DIM];
	int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);

	if(ndims != DIM) {
		cerr << " Wrong number of dimensions " << ndims << " - " << DIM << endl;
		exit(-22);
	}

	for(int i=0;i<3;i++){
		lbound[i]=0;
		ubound[i]=int(dims_out[DIM-(i+1)])-1;
	}
	data.resize(lbound,ubound);
  
	dataset.read(data, PredType::NATIVE_FLOAT );
	return true;
}



bool Hdf5iStream::Read3DMatrix(string ArrayName, NumMatrix<float,3> &data,
                               double *xb, double *dx)
{
	/* 
	   Routine to read 3D matrix from h5-file. Parameters are:
	   ArrayName -> Name of dataset w/o group name
	   data      -> NumMatrix-Array to hold data
	   xb        -> Array to hold values of lower bound positions
	   dx        -> Array to hold grid size
	*/

	const int DIM = 3;
	int lbound[3], ubound[3];
	float dummy[3];

	DataSet dataset = group->openDataSet( ArrayName );
	Attribute Origin = dataset.openAttribute("origin");
	Origin.read(PredType::NATIVE_FLOAT, dummy);
	xb[0] = double(dummy[0]);
	xb[1] = double(dummy[1]);
	xb[2] = double(dummy[2]);

	Attribute Delta = dataset.openAttribute("delta");
	Delta.read(PredType::NATIVE_FLOAT, dummy);
	dx[0] = double(dummy[0]);
	dx[1] = double(dummy[1]);
	dx[2] = double(dummy[2]);

	DataSpace dataspace = dataset.getSpace();

	int dimhdf = dataspace.getSimpleExtentNdims();
	if(dimhdf != DIM){
		cerr << " Wrong dimensionality of input data: " << endl;
		cerr << dimhdf << " " << DIM << endl;
		exit(-2);
	}
	hsize_t dims_out[DIM];
	dataspace.getSimpleExtentDims( dims_out, NULL);


	for(int i=0;i<3;i++){
		lbound[i]=0;
		ubound[i]=int(dims_out[DIM-(i+1)])-1;
	}
	data.resize(lbound,ubound);
  
	dataset.read(data, PredType::NATIVE_FLOAT );
	return true;
}



bool Hdf5Stream::Write3DMatrixSwap(string ArrayName, NumMatrix<double,3> &data,
                                   double *xb, double *dx)
{
	/* Routine to write NumMatix data in correct ordering to hdf5 file

	   Remarks: Only use if sufficient memory available - dummy
	   NumMatrix has to be created
	*/
	int hi[3];
	hi[0]=data.getHigh(2);
	hi[1]=data.getHigh(1);
	hi[2]=data.getHigh(0);
	int lo[3];
	lo[0] = data.getLow(2);
	lo[1] = data.getLow(1);
	lo[2] = data.getLow(0);

	NumMatrix<double,3> outdata(Index::set(lo[0],lo[1],lo[2]),
	                            Index::set(hi[0],hi[1],hi[2]));

	// Swapping data dimensions
	for(int i = lo[0]; i <= hi[0]; ++i){
		for(int j = lo[1]; j <= hi[1]; ++j){
			for(int k = lo[2]; k <= hi[2]; ++k){
				outdata(i,j,k) = data(k,j,i);
			}
		}
	}
  
	int mx[3];
	mx[2]= hi[0] - lo[0] + 1;
	mx[1]= hi[1] - lo[1] + 1;
	mx[0]= hi[2] - lo[2] + 1; 


	num+=1;
	if(num > NumEntries){
		return false;
	}
	FloatType datatype( PredType::NATIVE_DOUBLE );
	datatype.setOrder( H5T_ORDER_LE );
	const int DIM = 3;
	hsize_t DimsData[DIM];
	for(int q=0; q<DIM; ++q){
		DimsData[q]  = mx[q];
	}
	DataSpace dataspace( DIM, DimsData);


	// Supplying additional attributes for opendx input

	FloatType datatypefloat( PredType::NATIVE_DOUBLE );
	hsize_t DimsAttr = 3;
	DataSpace attrspace( 1, &DimsAttr );

	double Origin[3];
	double Delta[3];
	for(int q=0; q<3; ++q){
		Origin[q] = xb[q];
		Delta[q]  = dx[q];
	}


	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));

	// Writing attributes
	Attribute origin(dataset.createAttribute("origin", datatypefloat, attrspace));
	Attribute delta(dataset.createAttribute("delta", datatypefloat, attrspace));
	origin.write( datatypefloat, Origin);
	delta.write( datatypefloat, Delta);

	dataset.write( outdata, datatype );
	return true;
}



bool Hdf5Stream::Write2DMatrix(string ArrayName, NumMatrix<float,2> &data,
                               double *xb, double *dx, bool with_opendxinfo)
{
	/* Routine to write NumMatrix data in wrong ordering to hdf5 file.

	   Remarks:
     
	   On reading the hdf data the swapped dimensions have to be taken
	   into account
    
	*/
	int mx[2];
	mx[0]=data.getHigh(1) - data.getLow(1) + 1;
	mx[1]=data.getHigh(0) - data.getLow(0) + 1;
  
	AddDatasetName(ArrayName);

	num+=1;

	const int DIM = 2;
	FloatType datatype( PredType::NATIVE_FLOAT );
	datatype.setOrder( H5T_ORDER_LE );
  
	hsize_t DimsData[DIM];
	for(int q=0; q<DIM; ++q){
		DimsData[q]  = mx[q];
	}
	DataSpace dataspace( DIM, DimsData);
 
	// Supplying additional attributes for opendx input
	hsize_t DimsAttr = 2;
	DataSpace attrspace( 1, &DimsAttr );
	float Origin[DIM];
	float Delta[DIM];
	if( with_opendxinfo ) {
		for(int q=0; q<DIM; ++q){
			Origin[q] = float(xb[q]);
			Delta[q]  = float(dx[q]);
		}
	}
	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));

	// Writing attributes (if necessary)
	if( with_opendxinfo ) {
		Attribute origin(dataset.createAttribute("origin",datatype,attrspace));
		origin.write(datatype, Origin);

		Attribute delta(dataset.createAttribute("delta",datatype,attrspace));
		delta.write(datatype, Delta);
	}

	dataset.write(data, datatype);
	return true;
}



bool Hdf5Stream::WriteArray(string ArrayName, int *data, int max)
{
	num+=1;
	if(num > NumEntries){
		return false;
	}

	IntType datatype( PredType::NATIVE_INT );
	datatype.setOrder( H5T_ORDER_LE );

	const int DIM = 1;
	hsize_t DimsData[DIM];
	DimsData[0] = max;//sizeof(data)/sizeof(data[0]);

	DataSpace dataspace( DIM, DimsData);

	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));
	dataset.write( data, datatype );
	return true;
}



bool Hdf5Stream::WriteArray(int *data, int max)
{
	char ArrayName[256];
	sprintf(ArrayName,"data%3.3d",num+1);
	return WriteArray(ArrayName, data, max);
}





bool Hdf5Stream::WriteArray(string ArrayName, float *data, int max)
{
	num+=1;
	if(num > NumEntries){
		return false;
	}

	FloatType datatype( PredType::NATIVE_FLOAT );
	datatype.setOrder( H5T_ORDER_LE );

	const int DIM = 1;
	hsize_t DimsData[DIM];
	DimsData[0] = max;//sizeof(data)/sizeof(data[0]);

	DataSpace dataspace( DIM, DimsData);

	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));
	dataset.write( data, datatype);
	return true;
}



bool Hdf5Stream::WriteArray(float *data, int max)
{
	char ArrayName[256];
	sprintf(ArrayName,"data%3.3d",num+1);
	string AName = ArrayName;
	return WriteArray(AName, data, max);
}





bool Hdf5Stream::WriteArray(string ArrayName, double *data, int max)
{
	num+=1;
	if(num > NumEntries){
		return false;
	}

	FloatType datatype( PredType::NATIVE_DOUBLE );
	datatype.setOrder( H5T_ORDER_LE );

	const int DIM = 1;
	hsize_t DimsData[DIM];
	DimsData[0] = max;//sizeof(data)/sizeof(data[0]);

	DataSpace dataspace( DIM, DimsData);

	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));
	dataset.write( data, datatype );
	return true;
}




bool Hdf5Stream::WriteArray(double *data, int max)
{
	char ArrayName[256];
	sprintf(ArrayName,"data%3.3d",num+1);
	string AName = ArrayName;
	return WriteArray(AName, data, max);
}


bool Hdf5Stream::WriteNDArray(string ArrayName, float *data, int mx[], int dim)
{
	num+=1;
	if(num > NumEntries){
		return false;
	}

	FloatType datatype( PredType::NATIVE_FLOAT );
	datatype.setOrder( H5T_ORDER_LE );

	std::vector<hsize_t> DimsData(dim);
	for(int dir=0; dir<dim; ++dir) {
		DimsData[dir] = mx[dir];
	}

	DataSpace dataspace( 4, DimsData.data());

	DataSet dataset(group->createDataSet( ArrayName, datatype, dataspace ));
	dataset.write( data, datatype );
	return true;
}



bool Hdf5Stream::close()
{
	Attribute Info(group->openAttribute("Entries"));
	Info.write( PredType::NATIVE_INT, &num);

	delete group;
	delete hdf5file;
	this->open = false;
	return true;
}


bool Hdf5iStream::close()
{
	// delete group;
	// delete hdf5file;
	return true;
}


Hdf5Stream::~Hdf5Stream()
{
	
	if(open) {
		close();
	}
}


Hdf5iStream::~Hdf5iStream()
{
	delete group;
	delete hdf5file;
}


// Function: H5Object::doesAttrExist
///\brief test for existence of attribut
///\param name - IN: Name of the attribute
///\return true if attribute exists, false otherwise
///\exception none
// Programmer Kent Williams 2011 
bool Hdf5iStream::doesAttrExist( const char* name ) const {
	return( H5Aexists(group->getId(), name) > 0 ? true : false );
}
 
