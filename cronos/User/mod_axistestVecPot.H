#ifndef CRONOS_AXISTEST_VECPOT_H
#define CRONOS_AXISTEST_VECPOT_H = 1

#include "problem.H"
#include "randgen.H"
#include "vector.H"

class AxistestVecPot: public ProblemType {
public:
	AxistestVecPot (Data &);
	void init_fields (Data &, int [], int []);
private:
	short type_sfix;    // what type of singularity to fix
	double eps;         // something << dr_mean
	NumMatrix<int,1> mp;  // array of mapping partners
};


AxistestVecPot::AxistestVecPot(Data &gdata) : ProblemType(gdata)
{
	name = " grid singularity testing problem (VecPot) ";
	cerr << " Note: GEOM = " << GEOM << endl;
	
	/*  check if the current block contains a coordinate singularity
	 *  and do preparations if needed,
	 *  then set type_sfix to 0 -> no fix needed
	 *                        1 -> fix axis
	 *                        2 -> fix origin (still ToDo!)
	 *                        3 -> both
	 */
}


void AxistestVecPot::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{
	//	gdata.om[4].rename("A_x");
	// gdata.om[5].rename("A_y");
	// gdata.om[6].rename("A_z");
	double Rmax = gdata.xe[1];
	for (int k = ibeg[2]; k <= iend[2]; ++k) {
		for (int j = ibeg[1]; j <= iend[1]; ++j) {
			for (int i = ibeg[0]; i <= iend[0]; ++i) {
				REAL nrad = gdata.getCen_x(i)/gdata.xe[0];
				// REAL Az = (nrad < 0.5) ? cos(2.*M_PI*nrad)*gdata.xe[0]/(2.*M_PI) : 0.;
				REAL Az(0.);
				// if(nrad <= 0.5) {
				// 	Az = - 0.5*(gdata.get_x(i) - sin(4.*M_PI*nrad)*gdata.xe[0]/(4.*M_PI));
				// } else {
				// 	Az = -0.25*gdata.xe[0];
				// }
				// if(nrad < 0.) Az = - 0.5*(gdata.get_x(i));
				// if(nrad <= 0.5) {
				// 	Az = -0.5*(gdata.get_x(i) + sin(2.*M_PI*nrad)*gdata.xe[0]/(2.*M_PI));
				// } else {
				// 	Az = -0.25*gdata.xe[0];
				// }
				if(nrad <= 0.5) {
					Az = gdata.xe[0]*0.5 - gdata.getCen_x(i);
				} else {
					Az = 0.;
				}
				//				if(nrad < 0.) Az = - 0.5*(gdata.get_x(i));
				gdata.om[q_rho](i,j,k) = 1.0;

					// if(j==5 && k==0) {
					// 	cout << i << " " << gdata.getCen_x(i) << " ";
					// 	cout << Az << endl;
					// }

#if (GEOM == 2)
				REAL phi0 = gdata.getCen_y(j);
				double vx = 0.5;
				double vy = 0.;//0.5;
				gdata.om[q_sx](i,j,k) =  vx*cos(phi0) + vy*sin(phi0);
				gdata.om[q_sy](i,j,k) = -vx*sin(phi0) + vy*cos(phi0);

				gdata.om[q_Bx](i,j,k) =  0.;
				gdata.om[q_By](i,j,k) =  0.;
				gdata.om[q_Bz](i,j,k) =  Az;

#endif
				/*
				  bool ghost = ( (i < 0 || i > gdata.mx[0]) ||
				  (j < 0 || j > gdata.mx[1]) ||
				  (k < 0 || k > gdata.mx[2]) );
				*/
				if (ENERGETICS == FULL) {
					gdata.om[7](i,j,k) = 1.;
				}
			}
		}
	}
	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");
	// exit(2);

}






#endif // [AXISTESTHYDRO_H]
