#ifndef CRONOS_SAVE_DATA_H
#define CRONOS_SAVE_DATA_H = 1

#include "data.H"
#include "constants.H"
#include "matrix.H"
#include "fields_local.H"

class Saves {
public:
	Saves(const Data &);
	~Saves();
	void save_vars(const Data &,fields_1D &,
	               const int &, const int &, const int &);
	void retrieve(const Data &, fields_2D &, const int &, const int &);
private:
	NumMatrix<REAL,3>* save_field_x;
	NumMatrix<REAL,3>* save_field_y;
	NumMatrix<REAL,3>* save_field_z;
	int q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz;
};

#endif
