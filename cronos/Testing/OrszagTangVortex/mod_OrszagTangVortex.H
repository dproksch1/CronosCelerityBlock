#ifndef CRONOS_ORSZAG_TANG_VORTEX_H
#define CRONOS_ORSZAG_TANG_VORTEX_H

#include "problem.H"

// using namespace std;

/*
  Alle Benutzerklassen erben einfach von der Oberklasse Problemtype -
  die kann ich noch in specific.C unterbringen.

  Hier können beliebige benötigte Subroutinen eingebunden werden -
  damit stören sie nicht im Hauptcode
 */

class OrszagTangVortex: public ProblemType {
public:
	OrszagTangVortex (Data &);
	virtual void init_fields(Data &, int [], int []);
private:
	REAL get_init_sx(REAL x, REAL y);
	REAL get_init_sy(REAL x, REAL y);
	REAL get_init_Az(REAL x, REAL y);
	REAL B0, p0;
};

OrszagTangVortex::OrszagTangVortex(Data &gdata)  : ProblemType(gdata)
{
	B0 = value((char*)"B0");
	p0 = value((char*)"p0");
}


REAL OrszagTangVortex::get_init_sx(REAL x, REAL y) {
	return -sin(2. * M_PI * y);
}

REAL OrszagTangVortex::get_init_sy(REAL x, REAL y) {
	return sin(2. * M_PI * x);
}

REAL OrszagTangVortex::get_init_Az(REAL x, REAL y) {
	return B0/(4.*M_PI) * cos(4.*M_PI*x) + B0/(2.*M_PI) * cos(2.*M_PI*y);
}



void OrszagTangVortex::init_fields(Data &gdata, int ibeg[3], int iend[3])
{	

	for(int iz=ibeg[2]; iz<=iend[2]; ++iz) {
		for(int iy=ibeg[1]; iy<=iend[1]; ++iy) {
			for(int ix=ibeg[0]; ix<=iend[0]; ++ix) {
				
				REAL x_pos = gdata.getCen_x(ix);
				REAL y_pos = gdata.getCen_y(iy);

				gdata.om[q_rho](ix,iy,iz) = rho0;
				gdata.om[q_sx](ix,iy,iz) = get_init_sx(x_pos, y_pos);
				gdata.om[q_sy](ix,iy,iz) = get_init_sy(x_pos, y_pos);
				gdata.om[q_sz](ix,iy,iz) = 0.;


				gdata.om[q_Bx](ix,iy,iz) = 0.;
				gdata.om[q_By](ix,iy,iz) = 0.;

				// shift positions
				x_pos = gdata.getEdgL_x(ix);
				y_pos = gdata.getEdgL_y(iy);
				gdata.om[q_Bz](ix,iy,iz) = get_init_Az(x_pos, y_pos);


				gdata.om[q_Eges](ix,iy,iz) = p0 / (gamma - 1.);

			}
		}
	}



	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");


}





#endif
