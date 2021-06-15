#ifndef CRONOS_MRI_CYL_DISC_H
#define CRONOS_MRI_CYL_DISC_H = 1

#include "problem.H"

// using namespace std;

/*
  Alle Benutzerklassen erben einfach von der Oberklasse Problemtype -
  die kann ich noch in specific.C unterbringen.

  Hier können beliebige benötigte Subroutinen eingebunden werden -
  damit stören sie nicht im Hauptcode
*/

class MRICylDisc: public ProblemType {
public:
	MRICylDisc(Data &);
	virtual void init_fields(Data &, int [], int []);
	virtual void bc_User(Data &, NumMatrix<REAL,3> &, int, int, int, int);
	virtual void src_User(Data &, NumMatrix<REAL,3> [], NumMatrix<REAL,3> []);
	virtual void computeFluct(Data &, double &, double &);
private:
	REAL RandAmp, B0;
	RandomNumberGenerator* random;
};


MRICylDisc::MRICylDisc(Data &gdata) : ProblemType(gdata)
{
	name = " Cyldindrical Disc MRI ";
	// const string locname(" Cyldindrical Disc MRI ");
	// name = locname;
#if (GEOM != 2)
	cerr << " Error: Supposed to be a cylindircal problem: " << GEOM;
	cerr << endl;
	exit(2);
#endif
	this->RandAmp = value((char*)"Fluctuation_Amplitude");
// 	this->B0 = value((char*)"Initial_Magnetic_Field");
	this->B0 = 0.05513/8.0;

	REAL Lambda_crit = value((char*)"Critical_Wavelength");
	REAL R_crit = value((char*)"Critical_Radius");
	this->B0 = sqrt(rho0)*sqrt(15./16.)*Lambda_crit*pow(R_crit,-1.5)/(2.*pi);
// 	cout << " Using: " << B0*8 << endl;
// 	exit(2);

	long int seed = 42*gdata.rank;
	random = new RandomNumberGenerator(seed);
	for(int i=0; i<2000; ++i) {
		random->getRand();
	}
	// // Include normalisation:
	// TrafoNorm = new NormSound;
	// TrafoNorm->set_norms();

	//	cout << TrafoNorm << endl;

	// exit(2);

}



void MRICylDisc::src_User(Data &gdata, NumMatrix<REAL,3> nom[N_OMINT],
                          NumMatrix<REAL,3> nom_user[N_OMINT_USER]) {

	// Include gravitational force of central mass:
  
	for(int k=0; k<=gdata.mx[2]; ++k) {
		for(int j=0; j<=gdata.mx[1]; ++j) {
			for(int i=0; i<=gdata.mx[0]; ++i) {
				REAL r_cyl = gdata.getCen_x(i);
				REAL F_rad = 1./(sqr(r_cyl));

				nom[q_sx](i,j,k) += gdata.om[q_rho](i,j,k)*F_rad; 
	
				if(ENERGETICS == FULL) {
					if(gdata.om[q_sx].getName() == "v_x" &&
					   gdata.om[q_sy].getName() == "v_y") {
						nom[q_Eges](i,j,k) += gdata.om[q_rho](i,j,k)*gdata.om[q_sx](i,j,k)*F_rad;
					} else {
						cerr << " Error - has to be velocity " << endl;
						exit(2);
					}
				}
			}
		}
	}


}



void MRICylDisc::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			REAL phi = gdata.getCen_y(j);
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
				REAL r_cyl = gdata.getCen_x(i);

				// Density
				gdata.om[q_rho](i,j,k) = 1.;

				// Velocity
				gdata.om[q_sx](i,j,k) = 0.; // v_r
				gdata.om[q_sy](i,j,k) = sqrt(1./r_cyl); // v_phi
				gdata.om[q_sz](i,j,k) = 0.; // v_z

				// Random velocity fluctuations
				if(r_cyl > 1.5 && r_cyl < 3.5) {
					gdata.om[q_sx](i,j,k) += random->getRand(2*RandAmp) - RandAmp;
					gdata.om[q_sz](i,j,k) += random->getRand(2*RandAmp) - RandAmp;
				}
	
				// Vector potential:
				gdata.om[q_Bx](i,j,k) = 0.;
				gdata.om[q_By](i,j,k) = 0.;
				gdata.om[q_Bz](i,j,k) = 0.;

				// if(r_cyl > 2. && r_cyl < 3.) {
				// 	gdata.om[q_Bx](i,j,k) = -r_cyl*phi*B0;
				// }

				if(r_cyl > 2. && r_cyl < 3.) {
					// gdata.om[q_Bx](i,j,k) = -r_cyl*phi*B0;
					gdata.om[q_Bz](i,j,k) = B0;
				}


				// Thermal pressure:
				if(ENERGETICS == FULL) {
					gdata.om[q_Eges](i,j,k) = sqr(0.1)*gdata.om[q_rho](i,j,k)/gamma;
				}
			}
		}
	}

	// gdata.om[q_Bx].rename("A_x");
	// gdata.om[q_By].rename("A_y");
	// gdata.om[q_Bz].rename("A_z");

}



void MRICylDisc::bc_User(Data &gdata, NumMatrix<REAL,3> &omb,
                         int dir, int top, int q, int rim) 
{
  /*
    Include correct extrapolation for density and velocity for the
    initial conditions:
   */
  if(dir == 0) {
    int kmin = -1;
    if(q == 4) {
      kmin = -2;
    }
    if(top == 0) {
      for (int k = -rim; k <= gdata.mx[2]+rim; ++k){
	      for (int j = -rim; j <= gdata.mx[1]+rim; j++){
		      for (int i = kmin; i >= -rim; --i){
			      if(q == 2) {
				      REAL r_cyl0 = gdata.getCen_x(i+1);
				      REAL v_phi0 = gdata.om[q_sy](i+1,j,k)*sqrt(r_cyl0);
				      REAL r_cyl = gdata.getCen_x(i);
				      omb(i,j,k) = v_phi0/sqrt(r_cyl);
			      } else {
				      omb(i,j,k) = omb(i+1,j,k);
			      }
		      }
	      }
      }
    } else {
      for (int k = -rim; k <= gdata.mx[2]+rim; ++k){
	      for (int j = -rim; j <=gdata. mx[1]+rim; ++j){
		      for (int i = 1; i <= rim; ++i){
			      if (q == 2) {
				      REAL r_cyl0 = gdata.getCen_x(gdata.mx[0]+i-1);
				      REAL v_phi0 = gdata.om[q_sy](gdata.mx[0]+i-1,j,k)*sqrt(r_cyl0);
				      REAL r_cyl  = gdata.getCen_x(gdata.mx[0]+i);
				      omb(gdata.mx[0]+i,j,k) = v_phi0*sqrt(1./r_cyl);
			      } else {
				      omb(gdata.mx[0]+i,j,k) = omb(gdata.mx[0]+i-1,j,k);
			      }
		      }
	      }
      }
    }

  } else {
    cerr << " Only applicable for x-direction " << endl;
    exit(414);
  }
}




void MRICylDisc::computeFluct(Data &gdata, double &ekfluc, double &ebfluc)
{

	REAL Vol(0.);

	ekfluc = 0.;
	ebfluc = 0.;
	for (int k = 0; k < gdata.mx[2]; k++){
		for (int j = 0; j < gdata.mx[1]; j++){
			for (int i = 0; i < gdata.mx[0]; i++){

				REAL dV = gdata.get_CellVolume(i,j,k);
				
				double r_cyl = gdata.getCen_x(i);
				
				ekfluc += (sqr(gdata.om[q_sx](i,j,k)) +
				           sqr(gdata.om[q_sy](i,j,k) - sqrt(1./r_cyl)) +
				           sqr(gdata.om[q_sz](i,j,k)))*0.5*dV*gdata.om[q_rho](i,j,k);


				if(r_cyl > 2. && r_cyl < 3.) {
					ebfluc += (sqr(0.5*(gdata.om[q_Bx](i,j,k) + 
					                    gdata.om[q_Bx](i-1,j,k))) +
					           sqr(0.5*(gdata.om[q_By](i,j,k) + 
					                    gdata.om[q_By](i,j-1,k))) +
					           sqr(0.5*(gdata.om[q_Bz](i,j,k  ) - B0 +
					                    gdata.om[q_Bz](i,j,k-1) - B0)))*0.5*dV;
				} else {
					ebfluc += (sqr(0.5*(gdata.om[q_Bx](i,j,k) + 
					                    gdata.om[q_Bx](i-1,j,k))) +
					           sqr(0.5*(gdata.om[q_By](i,j,k) + 
					                    gdata.om[q_By](i,j-1,k))) +
					           sqr(0.5*(gdata.om[q_Bz](i,j,k  ) +
					                    gdata.om[q_Bz](i,j,k-1))))*0.5*dV;
				}


				Vol += dV;
			}
		}
	}

	ekfluc /= Vol;
	ebfluc /= Vol;

}

#endif
