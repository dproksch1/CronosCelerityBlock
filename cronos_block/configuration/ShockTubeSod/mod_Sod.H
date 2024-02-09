#ifndef CRONOS_MOD_DOD_H
#define CRONOS_MOD_SOD_H = 1

#include "problem.H"
#include "data.H"
#include "gridgen.H"
#include "utils.H"


class ShockTubeTestHydro: public ProblemType {
public:
	ShockTubeTestHydro(const Data &);
	virtual void init_fields(Data &, int [], int []);
	void bc_User(Data &, NumMatrix<double,3> &, int, int, int, int);
	void bc_User(Queue &, Data &, NumMatrix<double,3> &, int, int, int, int);
	void WriteToH5(Hdf5Stream &);
private:
	double rhoLInit, rhoRInit;
	double v0LInit, v0RInit;
	double VLXInit, VLYInit, VLZInit;
	double VRXInit, VRYInit, VRZInit;
	double PLInit, PRInit;
	double shockpos;
	int shockdir;
	bool isothermal;
};


ShockTubeTestHydro::ShockTubeTestHydro(const Data &gdata) : ProblemType(gdata) {

#if (FLUID_TYPE != CRONOS_HYDRO)
	cerr << endl;
	cerr << "!---------------------------------------------!"<< endl;
	cerr << "! ERROR: Need to use FLUID_TYPE CRONOS_HYDRO  !"<< endl;
	cerr << "!    --> correct in constants.H! Exiting...   !"<< endl;
	cerr << "!---------------------------------------------!"<< endl;
	cerr << endl;
	exit(2);
#endif

	if(ENERGETICS == FULL) {
		name = " Shock Tube problem (HD) 1D ";
		isothermal = false;
	} else {
		name = " Shock Tube problem (HD) 1D -- isothermal ";
		if(gamma > 1.) {
			cerr << " Has to be isothermal ";
			exit(2);
		}
		isothermal = true;
	}

	rhoLInit = value((char*)"rho0L");
	rhoRInit = value((char*)"rho0R");

	v0LInit = value((char*)"v0L");
	v0RInit = value((char*)"v0R");

	PLInit  = value((char*)"PL0");
	PRInit  = value((char*)"PR0");

	shockpos = value((char*)"shockposition");

	shockdir = static_cast<int>(value((char*)"shockdir"));

	if(shockdir==0) {
		VLXInit = v0LInit;
		VRXInit = v0RInit;
		VLYInit = 0.;
		VRYInit = 0.;
		VLZInit = 0.;
		VRZInit = 0.;
	} else if (shockdir==1) {
		VLXInit = 0.;
		VRXInit = 0.;
		VLYInit = v0LInit;
		VRYInit = v0RInit;
		VLZInit = 0.;
		VRZInit = 0.;
	} else {
		VLXInit = 0.;
		VRXInit = 0.;
		VLYInit = 0.;
		VRYInit = 0.;
		VLZInit = v0LInit;
		VRZInit = v0RInit;
	}

}

void ShockTubeTestHydro::init_fields(Data &gdata, int ibeg[3], int iend[3])
{
	double xPar;
	// gdata.om[q_Eges].rename("Temp");

	for (int k = ibeg[2]; k <= iend[2]; k++){
		double zz = gdata.getCen_z(k);
		for (int j = ibeg[1]; j <= iend[1]; j++){
			double yy = gdata.getCen_y(j);
			for (int i = ibeg[0]; i <= iend[0]; i++){
				double xx = gdata.getCen_x(i);

				if(shockdir == 0){
					xPar = xx;
				} else if (shockdir == 1){
					xPar= yy;
				} else {
					xPar = zz;
				}

				if(xPar < shockpos){
					gdata.om[q_rho](i,j,k) = rhoLInit;
					gdata.om[q_sx](i,j,k) = VLXInit;
					gdata.om[q_sy](i,j,k) = VLYInit;
					gdata.om[q_sz](i,j,k) = VLZInit;

					if(!isothermal) {
						gdata.om[q_Eges](i,j,k) = PLInit;
					}
				} else {
					gdata.om[q_rho](i,j,k) = rhoRInit;
					gdata.om[q_sx](i,j,k) = VRXInit;
					gdata.om[q_sy](i,j,k) = VRYInit;
					gdata.om[q_sz](i,j,k) = VRZInit;

					if(!isothermal) {
						gdata.om[q_Eges](i,j,k) = PRInit;
					}
				}
				if(!isothermal) {
					if(gdata.om[q_Eges].getName() == "Temp") {
						gdata.om[q_Eges](i,j,k) /= gdata.om[q_rho](i,j,k);
					} else {
						gdata.om[q_Eges](i,j,k) *= 1./(gamma-1.);
					}
				}
//				cout << " Velo: " << xPar << " " << xx << " " << gdata.om[q_sx](i,j,k) << " " << v0LInit << endl;
#if (CRSWITCH_DUAL_ENERGY == CRONOS_ON)
				gdata.om[q_Eadd](i,j,k) = gdata.om[q_Eges](i,j,k)*pow(gdata.om[q_rho](i,j,k),1-gamma);
#endif
	
			}
		}
	}
//	exit(3);
}



void ShockTubeTestHydro::bc_User(Data &gdata, NumMatrix<double,3> &omb,
		int dir, int top, int q, int rim) {
	if(dir == 0) {
		if(top == 1) {

			for (int k = -rim; k <= gdata.mx[2]+rim; ++k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = gdata.mx[0]+1; i <= gdata.mx[0]+B; ++i){
						if (q == q_rho) {      // density
							omb(i,j,k) = rhoRInit;
						} else if (q == q_sx) {
							omb(i,j,k) =  VRXInit;
						} else if (q == q_sy) {
							omb(i,j,k) =  VRYInit;
						} else if (q == q_sz) {
							omb(i,j,k) =  VRZInit;
						} else if (q == q_Eges) {
							omb(i,j,k) = PRInit/rhoRInit/(gamma-1.);
//							if(j==0 && k==0 && i==gdata.mx[0]+1)
//							cout << " Wert: " << PRInit << " " << omb(i,j,k) << " " << q << endl;
						} else {
							omb(i,j,k) = 0.;
						}
					}
				}
			}
		}
	}
}

void ShockTubeTestHydro::bc_User(Queue &queue, Data &gdata, NumMatrix<double,3> &omb,
		int dir, int top, int q, int rim) {

	if(dir == 0) {
		if(top == 1) {

			auto omRange = gdata.omSYCL[q].get_range();
			auto range = CelerityRange<3>(rim, omRange.get(1), omRange.get(2));

			double rhoRInit_local = this->rhoRInit;
			double VRXInit_local = this->VRXInit;
			double VRYInit_local = this->VRYInit;
			double VRZInit_local = this->VRZInit;
			double PRInit_local = this->PRInit;
			double gamma_local = gamma;
			double mx = gdata.mx[0];

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::all{}, celerity::write_only};

				cgh.parallel_for<class BCUserKernel>(range, [=](celerity::item<3> item) {

					size_t ix = item.get_id(0) + mx + rim + 1;
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					if (q == 0) {      // density
						om_acc[ix][iy][iz] = rhoRInit_local;
					} else if (q == 1) {
						om_acc[ix][iy][iz] =  VRXInit_local;
					} else if (q == 2) {
						om_acc[ix][iy][iz] =  VRYInit_local;
					} else if (q == 3) {
						om_acc[ix][iy][iz] =  VRZInit_local;
					} else if (q == 4) {
						om_acc[ix][iy][iz] = PRInit_local/rhoRInit_local/(gamma_local-1.);
					} else {
						om_acc[ix][iy][iz] = 0.;
					}
				});
			});
		}
	}
}

void ShockTubeTestHydro::WriteToH5(Hdf5Stream &h5out) {
	// Indicate that we are using a test case
	h5out.AddGlobalAttr("TestCase", "hydrodynamical shock tube test");

	// supply normalised values necessary for Parker solution
	// Use own group for that purpose
	hid_t group = h5out.AddGroup("Data/SodPars");

	// Supply left handed state variables
	h5out.AddGlobalAttr("rhoL", rhoLInit, group);
	h5out.AddGlobalAttr("uL", v0LInit, group);
	h5out.AddGlobalAttr("pL", PLInit, group);

	// Supply right handed state variables
	h5out.AddGlobalAttr("rhoR", rhoRInit, group);
	h5out.AddGlobalAttr("uR", v0RInit, group);
	h5out.AddGlobalAttr("pR", PRInit, group);


	// supply adiabatic exponent
	h5out.AddGlobalAttr("gamma", gamma, group);

	// Now close the new group
	h5out.CloseGroup(group);

}


#endif
