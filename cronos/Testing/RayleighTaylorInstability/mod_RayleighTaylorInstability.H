#ifndef CRONOS_RAYLEIGH_TAYLOR_INSTABILITY_H
#define CRONOS_RAYLEIGH_TAYLOR_INSTABILITY_H

#include "problem.H"
#include "randgen.H"

// using namespace std;

/*
  Alle Benutzerklassen erben einfach von der Oberklasse Problemtype -
  die kann ich noch in specific.C unterbringen.

  Hier können beliebige benötigte Subroutinen eingebunden werden -
  damit stören sie nicht im Hauptcode
 */

class RayleighTaylorInstability: public ProblemType {
public:
	RayleighTaylorInstability (Data &);
	virtual void init_fields(Data &, int [], int []);
	virtual void src_User(Data &, NumMatrix<REAL,3> [], NumMatrix<REAL,3> []);
	virtual void bc_User(Data &, NumMatrix<REAL,3> &, int, int, int, int);
	virtual void get_PhysFluxUser(Data &, cronos::vector<REAL> &, phys_fields_1D &, phys_fields_1D &, int, REAL);

private:
	REAL get_init(Data &, std::string, int, int, int);

	void add_stone_gardiner_perturbation(Data &, int, int, int);
	void add_mignone_perturbation(Data &, int, int, int);

	REAL rhoH, rhoL;
	REAL L, interface_pressure, adiabatic_exponent;
	REAL gx, gy, gz;
	#if (FLUID_TYPE == CRONOS_MHD)
	REAL Bx, By, Bz;
	#endif
	RandomNumberGenerator RandGen;
};

RayleighTaylorInstability::RayleighTaylorInstability(Data &gdata)  : ProblemType(gdata), RandGen(10)
{
	rhoH = value((char*)"rhoH");
	rhoL = value((char*)"rhoL");

	gx = value((char*)"gx");
	gy = value((char*)"gy");
	gz = value((char*)"gz");
	
	L = value((char*)"L");
	interface_pressure = value((char*)"interface_pressure");
	adiabatic_exponent = value((char*)"Adiabatic_exponent");
	
	

	#if (FLUID_TYPE == CRONOS_MHD)
	REAL g = sqrt(gx*gx + gy*gy + gz*gz);
	REAL Bc = sqrt((rhoH - rhoL) * g * L);

	Bx = value((char*)"Bx") * Bc;
	By = value((char*)"By") * Bc;
	Bz = value((char*)"Bz") * Bc;
	#endif
}




void RayleighTaylorInstability::init_fields(Data &gdata, int ibeg[3], int iend[3])
{	

	#if (FLUID_TYPE == CRONOS_MHD) 
		int qs[] = {q_rho, q_Eges, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz};
	#else
		int qs[] = {q_rho, q_Eges, q_sx, q_sy, q_sz};
	#endif


	for(int iz=ibeg[2]; iz<=iend[2]; ++iz) {
		for(int iy=ibeg[1]; iy<=iend[1]; ++iy) {
			for(int ix=ibeg[0]; ix<=iend[0]; ++ix) {

				
				for (int qi = 0; qi < gdata.fluid.get_N_OMINT(); ++qi) {
					gdata.om[qs[qi]](ix,iy,iz) = get_init(
						gdata, gdata.om[qs[qi]].getName(), ix, iy, iz);
				}

				for (int qi = 0; qi < gdata.fluid.get_N_OMINT_USER(); ++qi) {
					gdata.om_user[qi](ix,iy,iz) = get_init(
						gdata, gdata.om_user[qi].getName(), ix, iy, iz);
				}


				if (value((char*)"perturbation_type") == 0)
				{
				}
				else if (value((char*)"perturbation_type") == 1)
				{
					add_stone_gardiner_perturbation(gdata, ix, iy, iz);
				}
				else if (value((char*)"perturbation_type") == 2)
				{
					add_mignone_perturbation(gdata, ix, iy, iz);
				}
				else 
				{ 
					// throw std::invalid_argument("Unknown Perturbation Type");
				}



			}
		}
	}

}



void RayleighTaylorInstability::src_User(Data &gdata, NumMatrix<REAL,3> nom[N_OMINT],
                         NumMatrix<REAL,3> nom_user[N_OMINT_USER])
{
	for(int iz=0; iz<=gdata.mx[2]; ++iz) {
		for(int iy=0; iy<=gdata.mx[1]; ++iy) {
			for(int ix=0; ix<=gdata.mx[0]; ++ix) {
			
				// Adding linearized Gravity to the momentum equation
				nom[q_sx](ix,iy,iz) -= gdata.om[q_rho](ix,iy,iz) * gx;
				nom[q_sy](ix,iy,iz) -= gdata.om[q_rho](ix,iy,iz) * gy;
				nom[q_sz](ix,iy,iz) -= gdata.om[q_rho](ix,iy,iz) * gz;
				
				// Adding linearized Gravity to the energy equation
				nom[q_Eges](ix,iy,iz) -= gdata.om[q_rho](ix,iy,iz) * gx * gdata.om[q_sx](ix,iy,iz);
				nom[q_Eges](ix,iy,iz) -= gdata.om[q_rho](ix,iy,iz) * gy * gdata.om[q_sy](ix,iy,iz);
				nom[q_Eges](ix,iy,iz) -= gdata.om[q_rho](ix,iy,iz) * gz * gdata.om[q_sz](ix,iy,iz);
				
			}
		}
	}
}


void RayleighTaylorInstability::bc_User(Data &gdata, NumMatrix<REAL,3> &omb,
		int dir, int top, int q, int rim)
{
	/*
    	Include fixed boundary conditions for the z-Axis:
	Ghost cells will be filled with the initial values
	 */

	if(dir == 2) {

		if (top == 0) {
			for (int k = -rim; k < 0; ++k) {
				for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
					for (int i = -rim; i <= gdata.mx[0]+rim; ++i) {
						omb(i,j,k) = get_init(gdata, omb.getName(), i, j, k);
					}
				}
			}

		} else {
			for (int k = gdata.mx[2]+1; k <= gdata.mx[2]+rim; ++k) {
				for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
					for (int i = -rim; i <= gdata.mx[0]+rim; ++i) {
						omb(i,j,k) = get_init(gdata, omb.getName(), i, j, k);
					}
				}
			}
		}

	} 
}




inline void RayleighTaylorInstability::add_stone_gardiner_perturbation(Data &gdata, int ix, int iy, int iz)
{	
	REAL z_rel = gdata.getCen_z(iz)/L;

	if (abs(z_rel) < 0.5)
	{
		REAL A0 = 0.005;
		REAL R = 2.*RandGen.getRand() - 1.;
		gdata.om[q_sz](ix,iy,iz) += A0 * R * (1. + cos(2*M_PI*z_rel));
	}
}

inline void RayleighTaylorInstability::add_mignone_perturbation(Data &gdata, int ix, int iy, int iz)
{
	REAL x = gdata.getCen_x(ix);
	REAL y = gdata.getCen_y(iy);
	REAL z = gdata.getCen_z(iz);
	gdata.om[q_sz](ix,iy,iz) -= exp(-25.*(x*x+y*y))/(10.* pow(cosh(10.*z), 2) );
}


inline REAL RayleighTaylorInstability::get_init(Data &gdata, std::string name, int ix, int iy, int iz)
{
	if (name == "rho") {
		if (gdata.getCen_z(iz) >= 0) {
			return rhoH;
		} else  {
			return rhoL;
		}

	} else if (name == "Etherm") {
		REAL E_ges = interface_pressure;
		E_ges += gdata.om[q_rho](ix,iy,iz) * gdata.getCen_z(iz) * gz;
		E_ges /= adiabatic_exponent - 1.;
		return E_ges;
	
	} else if (name == "v_x" || name == "v_y" || name == "v_z") {
		return 0.;

	#if (FLUID_TYPE == CRONOS_MHD)
	} else if (name == "B_x") {
		return Bx;

	} else if (name == "B_y") {
		return By;

	} else if (name == "B_z") {
		return Bz;
	#endif

	} else if  (name == "om_user0") {
		if (gdata.getCen_z(iz) >= 0) {
			return 1*rhoH;
		} else  {
			return 0;
		}

	} else if  (name == "om_user1") {
		if (gdata.getCen_z(iz) >= 0) {
			return 0;
		} else  {
			return 1*rhoL;
		}

	} else {
		return 0.;
	}
}

void RayleighTaylorInstability::get_PhysFluxUser(Data &gdata, cronos::vector<REAL> &iPos, phys_fields_1D &pf, 
				phys_fields_1D &pfUser, int dir, REAL shift) {
	/*
	pfUser.flux[0].clear();
	for (int i = -1; i <= gdata.mx[dir]+1; ++i)
		pfUser.flux[0](i) = pfUser.uPri[0](i)*pf.uPri[1+dir](i);

	pfUser.flux[1].clear();
	for (int i = -1; i <= gdata.mx[dir]+1; ++i)
		pfUser.flux[1](i) = pfUser.uPri[1](i)*pf.uPri[1+dir](i);
	*/

	pfUser.flux[0].clear();
	for (int i = -1; i <= gdata.mx[dir]+1; ++i)
		pfUser.flux[0](i) = pfUser.uPri[0](i)*pf.uPri[1+dir](i)*pf.uPri[0](i);

	pfUser.flux[1].clear();
	for (int i = -1; i <= gdata.mx[dir]+1; ++i)
		pfUser.flux[1](i) = pfUser.uPri[1](i)*pf.uPri[1+dir](i)*pf.uPri[0](i);
}


#endif
