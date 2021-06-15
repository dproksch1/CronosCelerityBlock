#ifndef CRONOS_SEDOV_3D_SPHERICAL_H
#define CRONOS_SEDOV_3D_SPHERICAL_H = 1

#include "problem.H"
#include "randgen.H"

// using namespace std;



class Sedov3DSpherical: public ProblemType {
public:
	Sedov3DSpherical(Data &);
	void init_fields(Data &, int [], int []);
	void WriteToH5(Hdf5Stream &h5out);
private:
	REAL injected_Energy;
	REAL Sedov_Radius, Sedov_Volume;
	RandomNumberGenerator* random;
};


Sedov3DSpherical::Sedov3DSpherical(Data &gdata) : ProblemType(gdata)
{
	name = " Spherical Sedov Explosion ";
#if (GEOM != 3)
	cerr << " Error: Supposed to be a Spherical problem ";
	cerr << endl;
	exit(2);
#endif
	this->injected_Energy = value((char*)"Injected_Energy");
	this->Sedov_Radius = 3.02*gdata.dx[0];

	long int seed = 42;
	random = new RandomNumberGenerator(seed);
	for(int i=0; i<1000; ++i) {
		random->getRand();
	}
}



void Sedov3DSpherical::init_fields(Data &gdata, int ibeg[3], int iend[3])
{

	REAL Sedov_Volume(0);
	REAL full_Volume(0);
	REAL theta = gdata.getCen_y(0);
	REAL phi = gdata.getCen_z(0);
	for(int i=ibeg[0]; i<=iend[0];++i) {
		REAL r_sph = gdata.getCen_x(i);
		if(r_sph < Sedov_Radius) {
			// Add cell:
			Sedov_Volume += sin(theta)*gdata.dx[1]*gdata.dx[2]*sqr(r_sph)*gdata.dx[0];
			full_Volume += 4.*M_PI*sqr(r_sph)*gdata.dx[0];
		}
	}

	REAL num_radCells(0);
	for(int i=0; i<=iend[0];++i) {
		REAL r_sph_r = gdata.getEdgL_x(i+1);
		if(r_sph_r < Sedov_Radius) {
			num_radCells++;
			cout << i << " " << r_sph_r << " " << Sedov_Radius << " " << num_radCells << endl;
		} else {
			break;
		}
	}

	// Get effective radius:
	REAL r_Sedov = gdata.getEdgL_x(num_radCells);
	cout << " Nomma radius " << r_Sedov << endl;
	REAL Sedov_Volume_sph = 4.*M_PI/3.*r_Sedov*r_Sedov*r_Sedov;

//	Sedov_Volume = gdata.dx[0]*Number_Cells*4.*M_PI;
	if(gdata.rank == 0) {
		cout << " Initial volume: " << Sedov_Volume << " " << full_Volume << " " << sin(theta) << " " << Sedov_Volume_sph <<endl;
	}


	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
				REAL r_sph = gdata.getCen_x(i);

				gdata.om[q_rho](i,j,k) = rho0;
				gdata.om[q_sx](i,j,k) = 0.;
				gdata.om[q_sy](i,j,k) = 0.;
				gdata.om[q_sz](i,j,k) = 0.;
				gdata.om[q_Eges](i,j,k) = 1.e-5/(gamma-1.);

				if(r_sph < Sedov_Radius) {
//					gdata.om[q_Eges](i,j,k) += injected_Energy/Sedov_Volume;
//					gdata.om[q_Eges](i,j,k) += injected_Energy/full_Volume;
					gdata.om[q_Eges](i,j,k) += injected_Energy/Sedov_Volume_sph;
				}
				
			}
		}
	}
//	exit(2);

}


void Sedov3DSpherical::WriteToH5(Hdf5Stream &h5out) {
  // Indicate that we are using a test case
  h5out.AddGlobalAttr("TestCase", "Sedov Explosion");

  // supply normalised values necessary for Parker solution
  // Use own group for that purpose
  hid_t group = h5out.AddGroup("Data/SedovPars");

  // supply dimensionality of problem
  h5out.AddGlobalAttr("Dimension", 3, group);

  // supply adiabatic exponent
  h5out.AddGlobalAttr("Gamma", gamma, group);

  // supply external density
  h5out.AddGlobalAttr("rho_a", rho0, group);

  // supply injected energy
  h5out.AddGlobalAttr("E_inj", injected_Energy, group);

  // Now close the new group
  h5out.CloseGroup(group);

}



#endif
