#ifndef CRONOS_ParkerSpherical1D_iso_H
#define CRONOS_ParkerSpherical1D_iso_H = 1

#include "problem.H"
#include "randgen.H"

// using namespace std;
using namespace CRONOS_CONSTANTS;


class NormSolarWind : public normalisation {
public:
	virtual void setup();
};


void NormSolarWind::setup(
		)
{
	REAL gamma = value("Adiabatic_exponent");

	set_definition(LEN, SolarRadius);
	set_definition(MASS, 0.5 * HydrogenMass);
	set_definition(NUM_DENS, 1.e14 / cube(Meter));
	set_definition(TEMP, value("Norm_Temp") * Kelvin);

	// Velocity - here: speed of sound
	set_definition(VEL, sqrt(gamma * BoltzmannConstant * normQuantity(TEMP) / normQuantity(MASS)));
	set_definition(MAG_IND, sqrt(VacuumPermeability * normQuantity(NUM_DENS) * normQuantity(MASS)) * normQuantity(VEL));

	finish_setup();
}



class ParkerSph1DIso: public ProblemType {
public:
	ParkerSph1DIso(Data &);
	void init_fields(Data &, int [], int []);
	void bc_User(Data &, NumMatrix<REAL,3> &, int, int, int, int);
	void src_User(Data &, NumMatrix<REAL,3> [], NumMatrix<REAL,3> []);
	void WriteToH5(Hdf5Stream &);
private:
	REAL vr_Parker_isothSW(REAL);
	Quantity MStar, RStar, MassLossRate;
	Quantity RCrit, vCrit, rho0Crit;
	Quantity GTimesM;

	REAL MStar_num, RStar_num, MassLossRate_num;
	REAL RCrit_num, vCrit_num, rho0Crit_num;
	REAL GTimesM_num;

	//NumMatrix<double,3> bc_rho;
	//NumMatrix<double,3> bc_s_r;
	RandomNumberGenerator* random;
	REAL pos_comp1[DIM];
};


ParkerSph1DIso::ParkerSph1DIso(Data &gdata) : ProblemType(gdata)
{
	name = " isothermal 3D Parker wind ";
#if (GEOM != 3)
	cerr << " Error: Supposed to be a spherical problem ";
	cerr << endl;
	exit(2);
#endif


#if (ENERGETICS != ISOTHERMAL)
	cerr << " Energetics need to be isothermal for this test ";
#endif
	//Component 1
	MStar = value((char*)"Stellar_Mass") * SolarMass;
	RStar = value((char*)"Stellar_Radius") * SolarRadius;
	MassLossRate = value((char*)"Masslossrate") * SolarMass/Year;

	// Position of star usually at centre
	pos_comp1[0] =  0.;
	pos_comp1[1] =  0.;
	pos_comp1[2] =  0.;

	REAL gamma = value("Adiabatic_exponent");
	Quantity TempWind = value((char*)"Norm_Temp") * Kelvin;
	Quantity c_sound = sqrt(gamma * 2.*BoltzmannConstant*TempWind/HydrogenMass);

	// Critical velocity
	vCrit = sqrt(cs2) * c_sound;

	// Critical radius:
	RCrit = GravitationalConstant * MStar /(2. * sqr(vCrit));

	// Derive density from mass-loss rate:
	rho0Crit = MassLossRate / (4. * pi * sqr(RStar) * vCrit * RStar/RCrit);

	cout << " Mass-loss rate " <<   value((char*)"Masslossrate") << " MSun per year" << endl;
	cout << " Mass-loss rate " <<   MassLossRate << " kg/s " << endl;
	cout << " Density at stellar radius " << rho0Crit << endl;

	// Compute all necessary normalised variables
	TrafoNorm = new NormSolarWind();
    TrafoNorm->setup();

	MStar_num = TrafoNorm->get_num(MStar);
	RStar_num = TrafoNorm->get_num(RStar);
	MassLossRate_num = TrafoNorm->get_num(MassLossRate);

    RCrit_num = TrafoNorm->get_num(RCrit);
    vCrit_num = TrafoNorm->get_num(vCrit);
    rho0Crit_num = TrafoNorm->get_num(TrafoNorm->DENS, rho0Crit);

	// Normalised gravity
	GTimesM = GravitationalConstant * MStar;
	GTimesM_num = TrafoNorm->get_num(GTimesM);

	cout << " (normalised) stellar radius " << RStar << " (" << RStar_num << ")" << endl;
	cout << " (normalised) critical radius " << RCrit << " ("<<RCrit_num << ")" << endl;
	cout << " (normalised) critical velocity " << vCrit << " (" << vCrit_num << ")" << endl;
	cout << " (normalised) G*M " << GTimesM << " (" << GTimesM_num << ")" << endl;
	cout << " (normalised) mass loss rate " << MassLossRate << " (" << MassLossRate_num << ")" <<endl;
	cout << " (normalised) mass density at stellar radius " << rho0Crit << " (" << rho0Crit_num << ")" << endl;

	// Just a test for consistency:
	REAL rho0Crit_num2 = MassLossRate_num/(4. * pi * sqr(RStar_num) * vCrit_num * RStar_num / RCrit_num);
	cout << " dens via normalised quantities " << rho0Crit_num2 << endl;

//	// Compute density factor from mass-loss rate instead:


	// Density at solar radius (computed from the assumption of a
	// constant massloss rate)

	// Radius to which data are normalised
//	RNorm = value((char*)"Normalisation_Radius");

		
//	rho1 = MassLossRate/(4.*pi*sqr(RStar)*vr_Parker_isothSW(RNorm)*c_sound);

/*	bc_rho.resize(Index::set(-B,-B,-B),
	              Index::set(-1,gdata.mx[1]+B,gdata.mx[2]+B));
	bc_s_r.resize(Index::set(-B,-B,-B),
	              Index::set(-1,gdata.mx[1]+B,gdata.mx[2]+B));
*/
//	 exit(3);
}



void ParkerSph1DIso::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{	
	REAL v_2 = vr_Parker_isothSW(RStar_num);
//	REAL v_2 = vr_Parker_isothSW(0.1*RStarNorm);
//	REAL v_1 = vr_Parker_isothSW(RNorm);
//	REAL v_RNorm = vr_Parker_isothSW(RNorm*0.5);
	REAL v_RCrit_num = vr_Parker_isothSW(RCrit_num);
//	cout << v_RCrit_num << endl; exit(3);
	
	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {

			  REAL r_sph = gdata.getCen_x(i);

			  gdata.om[q_rho](i,j,k) = rho0Crit_num/cube(r_sph);
			  gdata.om[q_rho](i,j,k) = rho0Crit_num/sqr(r_sph);

			  gdata.om[q_sx](i,j,k) = v_RCrit_num*(r_sph)/(RCrit_num);///radial component;
			  if(r_sph > 2.*RCrit_num) {
				  gdata.om[q_sx](i,j,k) =2.*v_RCrit_num;

			  }

			  gdata.om[q_sx](i,j,k) = vr_Parker_isothSW(r_sph);


			  gdata.om[q_sy](i,j,k) = 0.;///r_sph1;

			  gdata.om[q_sz](i,j,k) = 0.;///r_sph1;
				
			}
		}
	}
}


void ParkerSph1DIso::src_User(Data &gdata, NumMatrix<REAL,3> nom[N_OMINT],
                         NumMatrix<REAL,3> nom_user[N_OMINT_USER]) {

	for(int k=0; k<=gdata.mx[2]; ++k) {
		for(int j=0; j<=gdata.mx[1]; ++j) {
			for(int i=0; i<=gdata.mx[0]; ++i) {
				REAL r_sph = gdata.getCen_x(i);
				
				REAL F_rad = GTimesM_num/(sqr(r_sph));
				
				if(i>-1 && i<4) {
					cout << " Force " << i << " " << F_rad << endl;
				}

				nom[q_sx](i,j,k) += gdata.om[q_rho](i,j,k)*F_rad;
				
			}
		}
	}
		

}

void ParkerSph1DIso::bc_User(Data &gdata, NumMatrix<REAL,3> &omb,
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

			static const int im[] = {    -1,   -1,   -1,   -1,   -2,   -1,   -1};
			static const int qi[] = {q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz};
			kmin = im[qi[q]];

			for (int k = -rim; k <= gdata.mx[2]+rim; ++k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = kmin; i >= -rim; --i){
						if(q != q_sx) {
							if (q == q_rho) {      // density
								REAL r_sph = gdata.getCen_x(i);
								omb(i,j,k) = std::abs(MassLossRate_num/(4.*M_PI*sqr(r_sph)*gdata.om[q_sx](i,j,k)));
								omb(i,j,k) = rho0Crit_num/sqr(r_sph);
//								omb(i,j,k) = rho0Crit_num;

								// Adjust the mass flux


							} else if (q == q_sy) {
								omb(i,j,k) = 0.;//omega*rad*sin(tet);
							} else {
								omb(i,j,k) = 0.;
							}
						} else {
							if(omb(i,j,k) < 1.e-2) {
								omb(i,j,k) = 1.e-2;
							} else {
								omb(i,j,k) = 2*omb(i+1,j,k) - omb(i+2,j,k);
							}
//							if(i<0) {
//								omb(i,j,k) = omb(0,j,k);
//							}
							REAL r_sph = gdata.getCen_x(i);
							omb(i,j,k) = vr_Parker_isothSW(r_sph);

////							double r_sph = gdata.getCen_x(i);
//							double r_sph_p2 = gdata.getCen_x(i+2);
//							double r_sph_p1 = gdata.getCen_x(i+1);
//
//							REAL rhoxv_p2 = gdata.om[q_rho](i+2,j,k)*gdata.om[q_sx](i+2,j,k)*sqr(r_sph_p2);
//							REAL rhoxv_p1 = gdata.om[q_rho](i+1,j,k)*gdata.om[q_sx](i+1,j,k)*sqr(r_sph_p1);
//							REAL rhoxv = 2*rhoxv_p1 - rhoxv_p2;
//
//
//
//							omb(i,j,k) = rhoxv/(gdata.om[q_rho](i,j,k)*sqr(r_sph));
//							if(omb(i,j,k) < 0.) omb(i,j,k) = 0.;
//							omb(i,j,k) = 0.;
						}
					}
				}
			}
		}
	} else if (dir == 1) {  // aligned field/flow at theta_{min|max}
		//                       (via reflecting bcons)
		//   A/S coefficents:      n  vr  vt  vp  Br  Bt  Bp   T

		const short as_coef[] = { +1, +1, -1, +1, +1, -1, +1, +1 };
		const short coe = as_coef[q];
		for (int k = -rim; k <= gdata.mx[2]+rim; k++) {
			for (int i = -rim; i <= gdata.mx[0]+rim; i++) {
				if (top == 0)        // theta_min
					if (q == q_By) {  //     B_theta
						omb(i,-1,k) = 0.;
						for (int dj = 2; dj <= rim; dj++)
							omb(i,-dj,k) = -1. * omb(i,dj-2,k);
					} else {          // not B_theta
						for (int dj = 1; dj <= rim; dj++)
							omb(i,-dj,k) = coe * omb(i,dj-1,k);
					}
				else                 // theta_max (top == 1)
					if (q == q_By) {  //     B_theta
						omb(i,gdata.mx[1],k) = 0.;
						for (int dj = 1; dj <= rim; dj++)
							omb(i,gdata.mx[1]+dj,k) = -1. * omb(i,gdata.mx[1]-dj,k);
					} else {          // not B_theta
						for (int dj = 1; dj <= rim; dj++)
							omb(i,gdata.mx[1]+dj,k) = coe * omb(i,gdata.mx[1]-dj+1,k);
					}
			}  // end i loop
		}    // end k loop
	}
}

REAL ParkerSph1DIso::vr_Parker_isothSW(REAL r_sph1)
{

	// Normalise distance from sun
	REAL r_norm = r_sph1/RCrit_num;

	/*
	  solve u2 - ln u2 = 4*(ln x + 1/x) - 3
	  
	  with u=v_r1/c_s  and x = r / r_crit, r_crit = 0.5*G*M/c_s2 = 0.5*gsurfacc*R2/c_s2
	  
	  using the Newton-Raphson method, i.e. we solve f(u)=0
	  with the formula unew=uold-f(uold)/f'(uold)
	  
	  We have 0=f(u)=4*(ln x + 1/x) -3 +ln u2 -u2
	  f'(u)=2/u-2u
	*/

	REAL rhs = 4.*(log(r_norm) + 1./r_norm) - 3.;

	REAL u=0.;


	/* Starting value for iterative Newton-Raphson method */
	
	if ( r_norm < 1.0) {
		u=0.000000000001; // Should be ok for most cases....
	}  else {
		u=1.1;
	}

	/* 50 iterations should be enough */

	for (int n=1; n<50; n++) {
		REAL unew = u - ( (log(sqr(u))-sqr(u)+rhs)/((2.0/u)-2.0*u) );
		u = unew;
	}

	// Multiply by normalised cs2 to find result.
	return u*sqrt(cs2);

} 

		
void ParkerSph1DIso::WriteToH5(Hdf5Stream &h5out) {
	// Indicate that we are using a test case
	h5out.AddGlobalAttr("TestCase", "isothermal Parker wind");

	// supply normalised values necessary for Parker solution
	// Use own group for that purpose
	hid_t group = h5out.AddGroup("Data/ParkerPars");

	// supply normalised critical radius:
	h5out.AddGlobalAttr("CriticalRadius", RCrit_num, group);

	// supply normalised speed of sound
	h5out.AddGlobalAttr("cs", sqrt(cs2), group);

	// Now close the new group
	h5out.CloseGroup(group);

}


#endif
