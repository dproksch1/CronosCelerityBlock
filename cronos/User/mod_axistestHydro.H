#ifndef CRONOS_AXISTEST_HYDRO_H
#define CRONOS_AXISTEST_HYDRO_H = 1

#include "problem.H"
#include "randgen.H"
#include "vector.H"

class AxistestHydro: public ProblemType {
public:
	AxistestHydro (Data &);
	void init_fields (Data &, int [], int []);
	void bc_User (Data &, NumMatrix<REAL,3> &, int, int, int, int);
	void src_User(Data &gdata,
	              NumMatrix<REAL,3> nom     [N_OMINT     ],
	              NumMatrix<REAL,3> nom_user[N_OMINT_USER]);
private:
	short type_sfix;    // what type of singularity to fix
	double eps;         // something << dr_mean
	NumMatrix<int,1> mp;  // array of mapping partners
};


AxistestHydro::AxistestHydro(Data &gdata) : ProblemType(gdata)
{
	name = " grid singularity testing problem (hydro) ";
	cerr << " Note: GEOM = " << GEOM << endl;
	
	/*  check if the current block contains a coordinate singularity
	 *  and do preparations if needed,
	 *  then set type_sfix to 0 -> no fix needed
	 *                        1 -> fix axis
	 *                        2 -> fix origin (still ToDo!)
	 *                        3 -> both
	 */
#if (GEOM != 1)  // nothing to do for Cartesian grids.
#ifdef parallel
	if (gdata.coords[1] != 0) {
		cerr << "error: Axis patch currently does not"
		     << "support >1 cores in phi dir. Aborting." << endl;
		exit(1);
	}
#endif // [parallel]
	this->eps = (1.e-5)*(gdata.xe-gdata.xb)/gdata.mx[1];
	this->mp.resize (Index::set(-B), Index::set(gdata.mx[1]+B));
	bool fix_axis = false;
#if (GEOM == 2)
	if (abs(gdata.getEdgL_x(0)) < eps) {  // only proceed if block touches r=0
		fix_axis = true;
		for (int i = -B; i < 0; ++i)  {// check for symmetric spacing in r
			if (abs(gdata.getCen_x(i)+gdata.getCen_x(-1-i)) > eps) {
				cerr << "cannot map: asymmetric radial spacing." << endl;
				exit(1);
			}
		}
		const double Pi2 = 2.*M_PI;
		for (int j = -B; j <= gdata.mx[1]+B; ++j) { // look for mapping p.s
			mp(j) = -100;   // default value to flag a failure
			double phi = gdata.getCen_y(j) + M_PI; // required mp position aim
			phi -= Pi2*int(phi/Pi2);    // keep phi in [0,2pi]
			if (gdata.xb[1] <= phi && phi <= gdata.xe[1]) { // mp inside grid
				for (int jm = 0; jm <= gdata.mx[1]; ++jm) { // loop candidates
					if (abs(gdata.getCen_y(jm) - phi) < eps) {   // found one!
						mp(j) = jm;
						break;
					}
				}
				if (mp(j) == -100) {  // mp needed, but none found
					cerr << "error: singularity present, "
					     << "but cells cannot be mapped." << endl;
					exit(1);
				}
			} else {
				mp(j) = -10;  // this means "no mp needed"
			}
		}
	}
	short fix_origin = false;   // (does not apply for GEOM=2)
	
#endif  // [GEOM == 2]
#if (GEOM == 3)
	cerr << "error: AxisPatch only for GEOM=2 so far." << endl;
	exit(1);
#endif
#if (GEOM == 2)
	this->type_sfix = (( (fix_axis  ) ? 1 : 0) +
	                   ( (fix_origin) ? 2 : 0));
#endif
#endif // [GEOM != 1]
}


void AxistestHydro::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{
	//	gdata.om[4].rename("A_x");
	// gdata.om[5].rename("A_y");
	// gdata.om[6].rename("A_z");
	double Rmax = gdata.xe[1];
	for (int k = ibeg[2]; k <= iend[2]; ++k) {
		for (int j = ibeg[1]; j <= iend[1]; ++j) {
			for (int i = ibeg[0]; i <= iend[0]; ++i) {
				gdata.om[q_rho](i,j,k) = 1.0;
#if (GEOM == 2)
				REAL nrad0 = gdata.getCen_x(i)/gdata.xe[0];
				REAL phi0 = gdata.getCen_y(j);
				REAL phiS = gdata.getEdgR_y(j);
				double vx = 0.5;
				double vy = 0.;//0.5;
				gdata.om[q_sx](i,j,k) =  vx*cos(phi0) + vy*sin(phi0);
				gdata.om[q_sy](i,j,k) = -vx*sin(phi0) + vy*cos(phi0);

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
	//	cout << " DFone " << endl;
}


void AxistestHydro::bc_User(Data &gdata, NumMatrix<REAL,3> &omb,
                    int dir, int top, int q, int rim)
{
	if (type_sfix == 0) return;
#if (GEOM == 2)
	if (dir == 0 && top == 0) {   // r_min
		/*  fill ghost cells with values from a computed grid cell for
		 *  cases in which both are located at the same physical position
		 *  (generic routine for both om and om_user fields)
		 *
		 *  NB: boundary cells without mapping partners are left unchanged!
		 */
		string qname = omb.getName(); 
		string subsp = qname.substr(qname.size()-2,2); // last two letters
		short si = (subsp == "_x" || subsp == "_y") ? -1 : +1; // sign inv.
		for (int k = -rim; k <= gdata.mx[2]+rim; ++k) {
			for (int j = -rim; j <= gdata.mx[1]+rim; ++j) {
				if (mp(j) != -10) {
					for (int i = -rim; i <= -1; i++) {
						omb(i,j,k) = si*omb(-1-i,mp(j),k);
					}
				} // mp(j)
			}
		}
		
	} // [if r_min]
#endif // [GEOM=2]
}


void AxistestHydro::src_User(Data &gdata,
                     NumMatrix<REAL,3> nom     [N_OMINT     ],
                     NumMatrix<REAL,3> nom_user[N_OMINT_USER])
{
// 	// --> make sure axis values are not evolved
// #if (GEOM == 2)
// 	if (type_sfix % 2 == 1)
// 		for (int k = 0; k <= gdata.mx[2]; k++)
// 			for (int j = 0; j <= gdata.mx[1]; j++)
// 				nom[q_Bx](-1,j,k) = 0.;   // remove NaN values
// 	// ToDo: vector potential case
// #endif
}

#endif // [AXISTESTHYDRO_H]
