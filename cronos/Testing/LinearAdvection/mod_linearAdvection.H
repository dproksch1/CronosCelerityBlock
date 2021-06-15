#ifndef CRONOS_ALFVEN_WINGS_SPHERICAL_H
#define CRONOS_ALFVEN_WINGS_SPHERICAL_H

#include "problem.H"

// using namespace std;

/*
  Alle Benutzerklassen erben einfach von der Oberklasse Problemtype -
  die kann ich noch in specific.C unterbringen.

  Hier können beliebige benötigte Subroutinen eingebunden werden -
  damit stören sie nicht im Hauptcode
 */

class LinearAdvection: public ProblemType {
public:
	LinearAdvection (Data &);
	virtual void init_fields(Data &, int [], int []);
	virtual void writePhystestInfo(Data &);
	virtual void src_User(Data &, NumMatrix<REAL,3> [], NumMatrix<REAL,3> []);
private:
	double get_FluctDataDirect(Data &);
	double set_disturbance(Data &, int ix, int iy, int iz);
	void WriteToH5(Hdf5Stream &h5out);
	REAL k_x, k_y;
	REAL rho_amp, v0;
	REAL velo_x, velo_y;
	REAL rho0;
	REAL deviation;
	REAL len0;
	int num_x, num_y;
};

LinearAdvection::LinearAdvection(Data &gdata)  : ProblemType(gdata)
{
#if (GEOM == CARTESIAN)
	name = " linear advect test ";
#else
	cerr << " Only Cartesian geometry is supported " << endl;
	exit(3);
#endif

	this->rho_amp = 0.05;
	this->rho_amp = value((char*)"Initial_Amplitude");
	num_x = value((char*)"modes_x");
	num_y = value((char*)"modes_y");

	this->v0 = 1.;


	// set wave number
	this->k_x = num_x*2.*M_PI/gdata.Lx[0];
	this->k_y = num_y*2.*M_PI/gdata.Lx[1];

	double norm = sqrt(sqr(num_x) + sqr(num_y));
	velo_x = v0*num_x/norm;
	velo_y = v0*num_y/norm;
  
	// Set correct size for array holding disturbance

	if(value_exists("BackgroundDensity")) {
		rho0 = value((char*)"BackgroundDensity");
	} else {
		rho0 = 1.;
	}

//	cF =
//	exit(3);

//	wave_type = 1; // shear Alfven waves
}





void LinearAdvection::init_fields(Data &gdata, int ibeg[3], int iend[3])
{

	for(int iz=ibeg[2]; iz<=iend[2]; ++iz) {
		for(int iy=ibeg[1]; iy<=iend[1]; ++iy) {
			for(int ix=ibeg[0]; ix<=iend[0]; ++ix) {

				//					gdata.om[q](ix,iy,iz) = 1.;

				gdata.om[q_rho](ix,iy,iz) = rho0 + set_disturbance(gdata,ix,iy,iz);
				gdata.om[q_sx](ix,iy,iz) = velo_x;
				gdata.om[q_sy](ix,iy,iz) = velo_y;
				gdata.om[q_sz](ix,iy,iz) = 0.;


			}
		}
	}

}


double LinearAdvection::set_disturbance(Data &gdata, int ix, int iy, int iz) {
	/*!
	 * Shear Alfven wave test taken from
	 * RyuEtAl1995
	 * the setup is given in Eq. (3.6) of that paper
	 * */
	REAL xPos = gdata.getCen_x(ix);
	REAL yPos = gdata.getCen_y(iy);


	REAL del_rho = rho_amp*sin(k_x*xPos + k_y*yPos);

	return del_rho;
}


double LinearAdvection::get_FluctDataDirect(Data &gdata) {
	double amp_vx(0.), amp_vy(0.), amp_vz(0.);
	double diff(0.);
	double num(0.);

	for( int ix=0; ix<=gdata.mx[0]; ++ix) {
		for( int iy=0; iy<=gdata.mx[1]; ++iy) {
			double val_init = rho0 + set_disturbance(gdata, ix, iy, 0);
			diff += sqr(gdata.om[q_sy](ix,iy,0) - val_init);
			num += 1.;
		}
	}

#ifdef parallel
	double numGlobal, diffGlobal;
	MPI_Allreduce(&num, &numGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Allreduce(&diff, &diffGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Barrier(gdata.comm3d);
	num = numGlobal;
	diff = diffGlobal;
#endif
	diff /= num;

//	cout << " amp: " << amplitude << endl;
//	exit(3);

	return diff;

}

void LinearAdvection::src_User(Data &gdata, NumMatrix<REAL,3> nom[N_OMINT],
                         NumMatrix<REAL,3> nom_user[N_OMINT_USER])
{
  for(int iz=0; iz<=gdata.mx[2]; ++iz) {
    for(int iy=0; iy<=gdata.mx[1]; ++iy) {
      for(int ix=0; ix<=gdata.mx[0]; ++ix) {

    	  // Using pure advection setup
    	  nom[q_sx](ix,iy,iz) = 0.;
    	  nom[q_sx](ix,iy,iz) = nom[q_rho](ix,iy,iz)*velo_x;
    	  nom[q_sy](ix,iy,iz) = nom[q_rho](ix,iy,iz)*velo_y;
    	  nom[q_sz](ix,iy,iz) = 0.;
//    	  gdata.om[q_sx](ix,iy,iz) = get_velo(gdata,ix);

      }
    }
  }
  //  get_rmsError(gdata);
}


void LinearAdvection::writePhystestInfo(Data &gdata) {
	// Here we use this to compute the fluctuation statistics at each relevant step


	double deviation = get_FluctDataDirect(gdata);

}


void LinearAdvection::WriteToH5(Hdf5Stream &h5out) {

	  // Indicate that we are using a test case
	  h5out.AddGlobalAttr("TestCase", "LinearAdvection");

	  // supply normalised values necessary for Parker solution
	  // Use own group for that purpose
	  hid_t group = h5out.AddGroup("Data/WavePars");

	  // suppy initial amplitude
	  h5out.AddGlobalAttr("Initial_Amplitude", rho_amp, group);

	  // Supply additional extent of domain
	  h5out.AddGlobalAttr("Len", len0, group);

	  // Background density
	  h5out.AddGlobalAttr("BackgroundDensity", rho0, group);

	  // Get deviation:
	  h5out.AddGlobalAttr("Error", deviation, group);

	  // Now close the new group
	  h5out.CloseGroup(group);

	}



#endif
