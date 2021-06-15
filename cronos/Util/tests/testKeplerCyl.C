#include "normalisation.H"
#include "KeplerOrbit.H"
#include <iostream>
#include <sstream>
#include <fstream>
#include <hdf5.h>

using namespace std;

// using namespace std;
using namespace CRONOS_CONSTANTS;


void WriteH5Dataset(hid_t group, string name, int num_steps, NumArray<double> &data) {
	// Append 1D array to h5 file
	hid_t datatype = H5Tcopy(H5T_NATIVE_DOUBLE);

	hsize_t DimsData = num_steps;
	hid_t dataspace = H5Screate_simple(1, &DimsData, NULL);

	hid_t dataset_id = H5Dcreate2(group, name.c_str(), datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Dwrite(dataset_id, datatype, H5S_ALL,
			H5S_ALL, H5P_DEFAULT, data);
	H5Dclose(dataset_id);

	H5Sclose(dataspace);
}


class NormStellarWind : public normalisation {
public:
  virtual void setup();
};


void NormStellarWind::setup(
    )
{
  REAL gamma = 5./3.;

  set_definition(LEN, SolarRadius);
  set_definition(MASS, 0.5 * HydrogenMass);
  set_definition(NUM_DENS, 1.e14 / cube(Meter));
  set_definition(TEMP, 10000 * Kelvin);

  // Velocity - here: speed of sound
  set_definition(VEL, sqrt(gamma * BoltzmannConstant * normQuantity(TEMP) / normQuantity(MASS)));
  set_definition(MAG_IND, sqrt(VacuumPermeability * normQuantity(NUM_DENS) * normQuantity(MASS)) * normQuantity(VEL));

  finish_setup();
}



int main() {
	ostringstream oss;
	oss << " Test, Test " << endl;
	cout << oss.str();

	// Get the norm

	NormStellarWind myNorm;
	myNorm.setup();
//	myNorm.setup();


//	const Quantity &massA, const Quantity &massB, const Quantity &semiMajor_Axis,
	Quantity massA = 120*SolarMass;
	Quantity massB = 30*SolarMass;
	Quantity period = 5.2*Year;

	double period_num = myNorm.get_num(period);
	oss.str("");
	oss << " timescale for full orbit: ";
	oss << period_num << endl;
	cout << oss.str();
	oss.str("");

	// Setting up a Keplerian orbit
	KeplerOrbit orbit = KeplerOrbit::newFromPeriod(massA, massB, period, 0.9, 0.5);

	KeplerOrbit *orbit2;
	orbit2 = new KeplerOrbit(KeplerOrbit::newFromPeriod(massA, massB, period, 0.9, 0.5));
	
	oss.clear();
	oss << orbit.get_semiMajorAxis(massA, massB, period) / AstronomicalUnit;
	oss << endl;
	cout << oss.str();


	NumArray<double> posA_num(3), posB_num(3);
	posA_num.clear();
	posB_num.clear();

	NumArray<double> veloA_num(3), veloB_num(3);
	veloA_num.clear();
	veloB_num.clear();


	hid_t hdf5file =  H5Fcreate("KeplerCyl.h5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	hid_t group    = H5Gcreate2(hdf5file, "/Orbit", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	int num_steps = 200;
	double del_t = period_num/(1.*num_steps);
	NumArray<double> h5_time(num_steps), h5_xA(num_steps), h5_yA(num_steps), h5_xB(num_steps), h5_yB(num_steps);
	NumArray<double> h5_vA_x(num_steps), h5_vA_y(num_steps), h5_vB_x(num_steps), h5_vB_y(num_steps);

	// Now, compute the position at a range of numerical times
	for(int i_time=0; i_time<num_steps; ++i_time) {
		ostringstream oss;
		double time_num = del_t*i_time;

		orbit.get_stellarPositionsCyl(myNorm, posA_num, posB_num, time_num);
		h5_time(i_time) = time_num;
		h5_xA(i_time) = posA_num(0);
		h5_yA(i_time) = posA_num(1);
		h5_xB(i_time) = posB_num(0);
		h5_yB(i_time) = posB_num(1);

		orbit.get_stellarVelocitiesCyl(myNorm, veloA_num, veloB_num, time_num);
		h5_vA_x(i_time) = veloA_num(0);
		h5_vA_y(i_time) = veloA_num(1);
		h5_vB_x(i_time) = veloB_num(0);
		h5_vB_y(i_time) = veloB_num(1);

		oss << " time " << time_num << " " << myNorm.get_phys(normalisation::TIME, time_num) << " ";
		oss << myNorm.get_phys(normalisation::TIME, time_num)/Year << " ";
		oss << endl;
		oss << "   ";
		oss << " A: (" << posA_num(0) << "," << posA_num(1) << ") ";
		oss << " B: (" << posB_num(0) << "," << posB_num(1) << ") ";
		oss << endl;
		cout << oss.str();

	}

	WriteH5Dataset(group, "time", num_steps, h5_time );
	WriteH5Dataset(group, "rhoA", num_steps, h5_xA );
	WriteH5Dataset(group, "phiA", num_steps, h5_yA );
	WriteH5Dataset(group, "rhoB", num_steps, h5_xB );
	WriteH5Dataset(group, "phiB", num_steps, h5_yB );
	WriteH5Dataset(group, "vAx", num_steps, h5_vA_x );
	WriteH5Dataset(group, "vAy", num_steps, h5_vA_y );
	WriteH5Dataset(group, "vBx", num_steps, h5_vB_x );
	WriteH5Dataset(group, "vBy", num_steps, h5_vB_y );


	H5Gclose(group);
	H5Fclose(hdf5file);

	return 0;
}
