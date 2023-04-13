#ifndef CRONOS_SEDOV_3D_H
#define CRONOS_SEDOV_3D_H = 1

#include "problem.H"
#include "randgen.H"

// using namespace std;



class Sedov3D: public ProblemType {
public:
	Sedov3D(Data &);
	void init_fields(Data &, int [], int []);
	void WriteToH5(Hdf5Stream &h5out);
private:
	double injected_Energy;
	double Sedov_Radius, Sedov_Volume;
	RandomNumberGenerator* random;
};


Sedov3D::Sedov3D(Data &gdata) : ProblemType(gdata)
{
	name = " Cartesian Sedov Explosion ";
#if (GEOM != 1)
	cerr << " Error: Supposed to be a Cartesian problem ";
	cerr << endl;
	exit(2);
#endif
	this->injected_Energy = value((char*)"Injected_Energy");
	this->Sedov_Radius = 3.5*gdata.dx[0];

	long int seed = 42;
	random = new RandomNumberGenerator(seed);
	for(int i=0; i<1000; ++i) {
		random->getRand();
	}
}



void Sedov3D::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{

	double Number_Cells(0);
	for(int k=ibeg[2]; k<=iend[2];++k) {
		double zval = k*gdata.dx[2];
		for(int j=ibeg[1]; j<=iend[1];++j) {
			double yval = j*gdata.dx[1];
			for(int i=ibeg[0]; i<=iend[0];++i) {
				double xval = i*gdata.dx[0];
				double rad = sqrt(sqr(xval) + sqr(yval) + sqr(zval));
				if(rad < Sedov_Radius) {
					// Add cell:
					Number_Cells += 1;
				}
			}
		}
	}
	Sedov_Volume = gdata.dx[0]*gdata.dx[1]*gdata.dx[2]*Number_Cells;
	if(gdata.rank == 0) {
		cout << " Initial volume: " << Sedov_Volume << endl;
	}


	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		double zz = gdata.getCen_z(k);
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			double yy = gdata.getCen_y(j);
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
				double xx = gdata.getCen_x(i);

				gdata.om[q_rho](i,j,k) = rho0;
				gdata.om[q_sx](i,j,k) = 0.;
				gdata.om[q_sy](i,j,k) = 0.;
				gdata.om[q_sz](i,j,k) = 0.;
				gdata.om[q_Eges](i,j,k) = 1.e-5/(gamma-1.);

				double rad = sqrt(sqr(xx) + sqr(yy) + sqr(zz));
				if(rad < Sedov_Radius) {
					gdata.om[q_Eges](i,j,k) += injected_Energy/Sedov_Volume;
				}
				
			}
		}
	}


}


void Sedov3D::WriteToH5(Hdf5Stream &h5out) {
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
