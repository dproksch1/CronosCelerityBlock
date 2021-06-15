#ifndef CRONOS_HYDRO_DISC_H
#define CRONOS_HYDRO_DISC_H = 1

#include "problem.H"
#include "randgen.H"

// using namespace std;



class HydroDisc: public ProblemType {
public:
	HydroDisc(Data &);
	void init_fields(Data &, int [], int []);
	void bc_User(Data &, NumMatrix<REAL,3> &, int, int, int, int);
	REAL c2_iso(Data &,const REAL &, const REAL &, const REAL &);
	void src_User(Data &, NumMatrix<REAL,3> [], NumMatrix<REAL,3> []);
	void get_PhysFluxUser(Data &gdata, cronos::vector<REAL> &iPos,
	                      phys_fields_1D &pf, phys_fields_1D &pfUser,
	                      int dir, REAL shift);
private:
	REAL Mass;
	RandomNumberGenerator* random;
};


HydroDisc::HydroDisc(Data &gdata) : ProblemType(gdata)
{
	name = " Hydrodynamical disc problem ";
#if (GEOM != 3)
	cerr << " Error: Supposed to be a spherical problem ";
	cerr << endl;
	exit(2);
#endif
	this->Mass = value((char*)"Central_Mass");

	long int seed = 42;
	random = new RandomNumberGenerator(seed);
	for(int i=0; i<1000; ++i) {
		random->getRand();
	}
}



void HydroDisc::src_User(Data &gdata, NumMatrix<REAL,3> nom[N_OMINT],
                         NumMatrix<REAL,3> nom_user[N_OMINT_USER]) {

  
	for(int k=0; k<=gdata.mx[2]; ++k) {
		for(int j=0; j<=gdata.mx[1]; ++j) {
			for(int i=0; i<=gdata.mx[0]; ++i) {
				REAL r_sph = gdata.getCen_x(i);
				REAL F_rad = Mass/(sqr(r_sph));

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



void HydroDisc::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			REAL theta = gdata.getCen_y(j);
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
				
				REAL r_sph = gdata.getCen_x(i);

				gdata.om[q_rho](i,j,k) = pow((r_sph*sin(theta)),-1.5)*exp((sin(theta) - 1.)/cs2);
				gdata.om[q_sx](i,j,k) = 0.;
				gdata.om[q_sy](i,j,k) = 0.;
				gdata.om[q_sx](i,j,k) += -0.0035 + random->getRand(0.007);
				gdata.om[q_sy](i,j,k) += -0.0035 + random->getRand(0.007);
	
				gdata.om[3](i,j,k) = sqrt((1./r_sph)*(1.-(2.5*cs2)/sin(theta)));

#if (OMS_USER == TRUE)
				if(i>gdata.mx[0]/4 && i <3*gdata.mx[0]/4) {
					if(j>gdata.mx[1]/4 && j <3*gdata.mx[1]/4) {
						if(k>gdata.mx[2]/4 && k <3*gdata.mx[2]/4) {
							gdata.om_user[0](i,j,k) = 1.;
						}
					}
				}
#endif
			}
		}
	}

#if (OMS_USER == TRUE)
	gdata.om_user[0].rename("Passive Scalar");
	set_q(-10, gdata.om_user[0]);
#endif
}


void HydroDisc::bc_User(Data &gdata, NumMatrix<REAL,3> &omb,
                        int dir, int top, int q, int rim) 
{
	/*
	  Include correct extrapolation for density and velocity for the
	  initial conditions:
	*/

	if(dir == 0) {
		int kmin = -1;
		if(q == q_Bx) {
			kmin = -2;
		}
		if(top == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; ++k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = kmin; i >= -rim; --i){
						if(q != 1 || omb(0,j,k) < 0) {
							omb(i,j,k) = omb(i+1,j,k);
						} else {
							omb(i,j,k) = 0;
						}
					}
				}
			}
		} else {
			for (int k = -rim; k <= gdata.mx[2]+rim; ++k){
				for (int j = -rim; j <=gdata. mx[1]+rim; ++j){
					for (int i = 1; i <= rim; ++i){
						if(q != 1 || omb(gdata.mx[0],j,k) > 0) {
							omb(gdata.mx[0]+i,j,k) = omb(gdata.mx[0]+i-1,j,k);
						} else {
							omb(i,j,k) = 0;
						}	  
					}
				}
			}
		}

	} else if (dir == 1) {

		if(top == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int j = -1; j >= -rim; --j){
						if(q != 2 || omb(i,0,k) < 0) {
							omb(i,j,k) = omb(i,j+1,k);
						} else {
							omb(i,j,k) = 0;
						}
					}
				}
			}
		} else {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int j = 1; j <= rim; j++){
						if(q != 2 || omb(i,gdata.mx[1],k) > 0) {
							omb(i,gdata.mx[1]+j,k) = omb(i,gdata.mx[1]+j-1,k);
						} else {
							omb(i,j,k) = 0;
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


void HydroDisc::get_PhysFluxUser(Data &gdata,
                                 cronos::vector<REAL> &iPos,
                                 phys_fields_1D &pf,
                                 phys_fields_1D &pfUser,
                                 int dir, REAL shift) {

	pfUser.flux[0].clear();
	
	// Example: transport of a passive tracer-field
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){
		
		pfUser.flux[0](i) = pfUser.uPri[0](i)*pf.uPri[1+dir](i);

	}

}
	

REAL HydroDisc::c2_iso(Data &gdata, 
                       const REAL &ii, const REAL &jj, const REAL &kk)
{
	REAL rad = gdata.getCen_x(ii);
	REAL theta = gdata.getCen_y(jj);
  
	REAL c2_loc = cs2/(rad*sin(theta)); // square of local sound speed
	return c2_loc;
}
		

#endif
