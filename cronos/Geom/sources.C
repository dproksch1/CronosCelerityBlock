#include <stdlib.h>
#include "sources.H"

SourceTerms::SourceTerms(const Data &gdata, const ProblemType &Problem,
                         bool IntegrateA, int thermal) {
	this->IntegrateA = IntegrateA;
	if(thermal==0) {
		this->thermal = false;
	} else {
		this->thermal = true;
	}

#if (GEOM == CARTESIAN) 
	geom = "Cartesian";
#elif (GEOM == CYLINDRICAL) 
	geom = "Cylindrical";
	facShift = gdata.dx[0]/12.;
#elif (GEOM == SPHERICAL)
	geom = "Spherical";
	facShift = 6./sqr(gdata.dx[0]);
#endif
	eos = new EquationOfState(Problem);

#ifdef GEOM_SOURCE_CORRECTION
	Limiter = new limiter;
#endif

#if (USE_COROTATION == CRONOS_ON)
	omegaZ = value((char*)"omegaZ");
#endif
}

SourceTerms::~SourceTerms() {
	delete eos;
}


// Nothing to be done for Cartesian geometry
#if (GEOM == CARTESIAN) 
void SourceTerms::src_Geom_Hydro(Data &gdata, ProblemType &Problem,
                                 NumMatrix<REAL,3> nom[], const CronosFluid &) {
	return;
}

#ifdef PHYSDISS
void SourceTerms::src_Geom_Viscosity(Data &gdata, ProblemType &Problem,
                                     NumMatrix<REAL,3> nom[], const CronosFluid &) {
	return;
}

void SourceTerms::mod_Geom_Visc_WE(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[], REAL tau[],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	return;
}

void SourceTerms::mod_Geom_Visc_SN(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[], REAL tau[],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	return;
}

void SourceTerms::mod_Geom_Visc_BT(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[], REAL tau[],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	return;
}

#endif
#endif



// Cylindrical geometry
#if (GEOM == CYLINDRICAL)

// Shift of collocation point for reconstruction polynomial
REAL SourceTerms::shift_Geom_WE(Data &gdata, 
                                int ii, int jj, int kk) {
#if (NON_LINEAR_GRID == CRONOS_ON)
	facShift = gdata.getCen_dx(0,ii)/12.;
#endif
	return facShift/gdata.getCen_x(ii);
}

REAL SourceTerms::shift_Geom_SN(Data &gdata, 
                                int ii, int jj, int kk) {
	return 0.;
}

REAL SourceTerms::shift_Geom_BT(Data &gdata, 
                                int ii, int jj, int kk) {
	return 0.;
}




void SourceTerms::src_Geom_Hydro(Data &gdata, ProblemType &pr,
                                 NumMatrix<REAL,3> nom[], const CronosFluid &fluid) {

	// Set indices:
	int q_rho = fluid.get_q_rho();
	int q_sx = fluid.get_q_sx();
	int q_sy = fluid.get_q_sy();
	int q_sz = fluid.get_q_sz();
	int q_Bx = fluid.get_q_Bx();
	int q_By = fluid.get_q_By();
	int q_Bz = fluid.get_q_Bz();
	int q_Eges = fluid.get_q_Eges();
	int q_Eadd = fluid.get_q_Eadd();

	int fluid_type = fluid.get_fluid_type();

	if(gdata.om[q_sx].getName() != "v_x" ||
	   gdata.om[q_sy].getName() != "v_y" ||
	   gdata.om[q_sz].getName() != "v_z") {
		throw CException(" Need characteristic variables for geom sources ");
	}

#ifdef GEOM_SOURCE_CORRECTION
	REAL fac = 1./12.;
#endif

	REAL B_r(0.), B_phi(0.), B_z(0.);
	for (int k = 0; k <= gdata.mx[2]; ++k){
		for (int j = 0; j <= gdata.mx[1]; ++j){
			for (int i = 0; i <= gdata.mx[0]; ++i){

				REAL r_cyl_inv = 1./(gdata.getCen_x(i));
				REAL rho(gdata.om[q_rho](i,j,k));
				REAL pres(0.);

				if(ENERGETICS == FULL) {
					if(thermal) {
						pres = (pr.gamma - 1.)*gdata.om[q_Eges](i,j,k);
					} else {
						pres = rho*gdata.om[q_Eges](i,j,k);
					}
				} else {
					pres = eos->pressure(gdata, pr, 
					                     gdata.om[q_rho](i,j,k), i, j, k);
				}


				REAL v_phi(gdata.om[q_sy](i,j,k));

#if (USE_COROTATION == CRONOS_ON)
				REAL r_cyl = gdata.getCen_x(i);
				REAL v_phiCorot = v_phi - r_cyl*omegaZ;
#endif

				//#if (FLUID_TYPE == CRONOS_MHD)
				if(fluid_type == CRONOS_MHD) {
				  B_r = 0.5*(gdata.om[q_Bx](i,j,k) +
							 gdata.om[q_Bx](i-1,j,k));
				  B_phi = 0.5*(gdata.om[q_By](i,j,k) +
							   gdata.om[q_By](i,j-1,k));
				  B_z = 0.5*(gdata.om[q_Bz](i,j,k) +
							 gdata.om[q_Bz](i,j,k-1));
				// } else {
				//   B_r = 0.;
				//   B_phi = 0.;
				//   B_z = 0.;
//					//#else
//						REAL B_r(0.), B_phi(0.), B_z(0.);
						//#endif
				}

				// if(i==gdata.mx[0]-10 && j==gdata.mx[1]/8 && k==0) {
				// 	REAL phi = gdata.getCen_y(j);
				// 	cout << " Sources: " << endl << "      ";
				// 	cout << B_r << " " << B_phi << " " << B_z << " ";
				// 	cout << sqrt(sqr(B_r) + sqr(B_phi)) << endl;
				// 	cout << " should: " << endl;
				// 	cout << "      ";
				// 	cout << cos(phi) << " " << -sin(phi) << endl;
				// 	cout << " more: " << endl;
				// 	cout << gdata.om[q_Bx](i,j,k) << " ";
				// 	cout << gdata.om[q_Bx](i-1,j,k) << endl;
				// 	cout << " alternative " << endl;
				// 	REAL phiL = gdata.getEdgL_y(j);
				// 	REAL phiR = gdata.getEdgR_y(j);
				// 	REAL fus = cos(phi)*0.5*(gdata.om[q_By](i,j,k)/cos(phiL) +
				// 	                         gdata.om[q_By](i,j-1,k)/cos(phiR));
				// 	cout << fus << endl;
				// 	B_phi = fus;
					
				// }

#if (USE_COROTATION == CRONOS_ON)
				nom[q_sx](i,j,k) -=  (rho*v_phiCorot*v_phi + pres +
						0.5*(sqr(B_r) + sqr(B_z) - sqr(B_phi)))*r_cyl_inv;
#else
				nom[q_sx](i,j,k) -=  (rho*sqr(v_phi) + pres +
						0.5*(sqr(B_r) + sqr(B_z) - sqr(B_phi)))*r_cyl_inv;
#endif


#ifdef GEOM_SOURCE_CORRECTION
				// Compute second order correction terms for the geometrical
				// source terms
				if(pr.mag) {
					if(gdata.rank == 0) {
						cerr << " Geometrical correction not implemented for mag field ";
					}
					exit(-3);
				}
#if (USE_COROTATION == CRONOS_ON)
					if(gdata.rank==0) {
						cerr << " Second order geom source correction not implemented for co-rotation " << endl;
					}
					exit(-3;)
#endif

				// Compute local gradients
				REAL drhodrp  =  gdata.om[q_rho](i+1,j  ,k  ) - gdata.om[q_rho](i  ,j  ,k  );
				REAL drhodr0  = (gdata.om[q_rho](i+1,j  ,k  ) -
				                 gdata.om[q_rho](i-1,j  ,k  ))*0.5;
				REAL drhodrm  =  gdata.om[q_rho](i  ,j  ,k  ) - gdata.om[q_rho](i-1,j  ,k  );
        
				REAL dvrdrp   =  gdata.om[q_sx](i+1,j  ,k  ) - gdata.om[q_sx](i  ,j  ,k  );
				REAL dvrdr0   = (gdata.om[q_sx](i+1,j  ,k  ) -
				                 gdata.om[q_sx](i-1,j  ,k  ))*0.5;
				REAL dvrdrm   =  gdata.om[q_sx](i  ,j  ,k  ) - gdata.om[q_sx](i-1,j  ,k  );
        
				REAL dvphidrp =  gdata.om[q_sy](i+1,j  ,k  ) - gdata.om[q_sy](i  ,j  ,k  );
				REAL dvphidr0 = (gdata.om[q_sy](i+1,j  ,k  ) -
				                 gdata.om[q_sy](i-1,j  ,k  ))*0.5;
				REAL dvphidrm =  gdata.om[q_sy](i  ,j  ,k  ) - gdata.om[q_sy](i-1,j  ,k  );
        
				// Get gradient using the limiter:
				REAL drhodr  = Limiter->compute(drhodrp,  drhodr0,  drhodrm);
				REAL dvrdr   = Limiter->compute(dvrdrp,   dvrdr0,   dvrdrm);
				REAL dvphidr = Limiter->compute(dvphidrp, dvphidr0, dvphidrm);

				nom[q_sx](i,j,k) -= -(rho*sqr(dvphidr) +
				                         2.*drhodr*v_phi*dvphidr)*fac*r_cyl_inv;

#endif

#if (USE_ANGULAR_MOMENTUM == FALSE)
				REAL v_r(gdata.om[q_sx](i,j,k));
				//				nom[q_sy](i,j,k) -= (-rho*v_r*v_phi + 0.5*(B_r*B_phi))*r_cyl_inv;
#if (USE_COROTATION == CRONOS_ON)
				nom[q_sy](i,j,k) -= (-rho*v_r*v_phiCorot + (B_r*B_phi))*r_cyl_inv;
#else
				nom[q_sy](i,j,k) -= (-rho*v_r*v_phi + (B_r*B_phi))*r_cyl_inv;
#endif

#ifdef GEOM_SOURCE_CORRECTION
				nom[q_sy](i,j,k) -= -(rho*dvrdr*dvphidr + drhodr*v_r*dvphidr +
				                         drhodr*dvrdr*v_phi)*fac*r_cyl_inv;
#endif
#endif

			}
		}
	}


}


#ifdef PHYSDISS
void SourceTerms::src_Geom_Viscosity(Data &gdata, ProblemType &pr,
                                     NumMatrix<REAL,3> nom[N_OMINT], const CronosFluid &fluid) {

	// Set indices:
	int q_rho = fluid.get_q_rho();
	int q_sx = fluid.get_q_sx();
	int q_sy = fluid.get_q_sy();
	int q_sz = fluid.get_q_sz();


	if(gdata.om[q_sx].getName() != "v_x" ||
	   gdata.om[q_sy].getName() != "v_y" ||
	   gdata.om[q_sz].getName() != "v_z") {
		throw CException(" Need characteristic variables for geom sources ");
	}

	Pot & rho   = gdata.om[q_rho];
	Pot & v_r   = gdata.om[q_sx];
	Pot & v_phi = gdata.om[q_sy];
	Pot & v_z   = gdata.om[q_sz];


	for (int k = 0; k <= gdata.mx[2]; ++k){
		for (int j = 0; j <= gdata.mx[1]; ++j){
			for (int i = 0; i <= gdata.mx[0]; ++i){
				REAL r_cyl_inv = 1./(gdata.getCen_x(i));

				double mu = pr.nu(gdata,i,j,k)*rho(i,j,k);

				//	REAL T_phiphi = 2*dv_phidphi - 2/3*div_v

				nom[q_sx](i,j,k) -= -mu*(2.*(v_phi(i,j+1,k) -
				                                v_phi(i,j-1,k))*gdata.getCen_hx(1,j)/gdata.getCen_h1(i,j,k) -
				                            ((v_r(i+1,j,k) - v_r(i-1,j,k))*gdata.getCen_hx(0,i) +
				                             (v_z(i,j,k+1) - v_z(i,j,k-1))*gdata.getCen_hx(1,j)) +
				                            (2.*v_r(i,j,k)*r_cyl_inv))*twothirds*r_cyl_inv;


				nom[q_sy](i,j,k) -= mu*((v_phi(i+1,j,k) -
				                            v_phi(i-1,j,k))*gdata.getCen_hx(0,i) -
				                           (v_phi(i,j,k)*r_cyl_inv) +
				                           (v_r(i,j+1,k) - 
				                            v_r(i,j-1,k))*gdata.getCen_hx(2,k)/gdata.getCen_h1(i,j,k))*r_cyl_inv;


			}
		}
	}

}


void SourceTerms::mod_Geom_Visc_WE(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[N2OMINT], REAL tau[N_OMINT],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	cerr << " sources::mod_Geom_Visc_WE - not adapted for multifluid yet " << endl;
	exit(3);
	REAL r_cyl_inv = 1./gdata.getCen_h1(i,j,k);

	// tau_rr
	tau[1] += -omWE[Problem.q_sx+dir]*r_cyl_inv*twothirds;
	// tau_rphi
	tau[2] += -omWE[Problem.q_sy+dir]*r_cyl_inv;

	return;
}


void SourceTerms::mod_Geom_Visc_SN(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[N2OMINT], REAL tau[N_OMINT],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	cerr << " sources::mod_Geom_Visc_SN - not adapted for multifluid yet " << endl;
	exit(3);
	REAL r_cyl_inv = 1./gdata.getCen_h1(i,j,k);

	// tau_phir
	tau[1] += -omWE[Problem.q_sy+dir]*r_cyl_inv;
	// tau_phiphi
	tau[2] +=  (2.*omWE[Problem.q_sx+dir]*r_cyl_inv)*twothirds;

	return;
}


void SourceTerms::mod_Geom_Visc_BT(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[N2OMINT], REAL tau[N_OMINT],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	cerr << " sources::mod_Geom_Visc_BT - not adapted for multifluid yet " << endl;
	exit(3);
	REAL r_cyl_inv = 1./gdata.getCen_h1(i,j,k);

	// tau_zz
	tau[3] += -omWE[Problem.q_sx+dir]*r_cyl_inv*twothirds;

	return;
}
#endif

#endif
#if(GEOM == SPHERICAL)

// Shift of collocation point for reconstruction polynomial
REAL SourceTerms::shift_Geom_WE(Data &gdata, int ii) {
#if (NON_LINEAR_GRID == CRONOS_ON)
	facShift = 6./sqr(gdata.getCen_dx(0,ii));
#endif
	return gdata.getCen_x(ii)/(0.5 + facShift*sqr(gdata.getCen_x(ii)));
}

REAL SourceTerms::shift_Geom_SN(Data &gdata, int jj) {
	return 0.;
}

REAL SourceTerms::shift_Geom_BT(Data &gdata, int kk) {
	return 0.;
}


void SourceTerms::src_Geom_Hydro(Data &gdata, ProblemType &pr,
                                 NumMatrix<REAL,3> nom[N_OMINT], const CronosFluid &fluid) {


	// Set indices:
	int q_rho = fluid.get_q_rho();
	int q_sx = fluid.get_q_sx();
	int q_sy = fluid.get_q_sy();
	int q_sz = fluid.get_q_sz();
	int q_Bx = fluid.get_q_Bx();
	int q_By = fluid.get_q_By();
	int q_Bz = fluid.get_q_Bz();
	int q_Eges = fluid.get_q_Eges();
	int q_Eadd = fluid.get_q_Eadd();

	int fluid_type = fluid.get_fluid_type();

	if(gdata.om[q_sx].getName() != "v_x" ||
	   gdata.om[q_sy].getName() != "v_y" ||
	   gdata.om[q_sz].getName() != "v_z") {
		throw CException(" Need characteristic variables for geom sources ");
	}

	REAL B_r(0.), B_theta(0.), B_phi(0.);
	for (int k = 0; k <= gdata.mx[2]; ++k){
		for (int j = 0; j <= gdata.mx[1]; ++j){
			REAL theta  = gdata.getCen_y(j);
			REAL thetap = gdata.getEdgR_y(j);
			REAL thetam = gdata.getEdgL_y(j);

			//      cot = cos(yy)/sin(yy);
			// Numerical approximation for cotanges of theta:
			REAL cot = 2.*(sin(thetap)-sin(thetam))/sin(theta)*gdata.getCen_hx(1,j);

			for (int i = 0; i <= gdata.mx[0]; ++i){
				REAL r_inv = 1./(gdata.getCen_x(i));
				REAL rho(gdata.om[q_rho](i,j,k));
				REAL pres(0.);

				if(ENERGETICS == FULL) {
					if(thermal) {
						pres = (pr.gamma - 1.)*gdata.om[q_Eges](i,j,k);
					} else {
						pres = rho*gdata.om[q_Eges](i,j,k);
					}
				} else {
					//          pres = pressure(gdata.om[q_rho](i,j,k));
					pres = eos->pressure(gdata, pr, gdata.om[q_rho](i,j,k),
					                     i, j, k);
				}


				REAL v_r(gdata.om[q_sx](i,j,k));
				REAL v_theta(gdata.om[q_sy](i,j,k));
				REAL v_phi(gdata.om[q_sz](i,j,k));

//#if (FLUID_TYPE == CRONOS_MHD)
				if(fluid_type == CRONOS_MHD) {
				  B_r = 0.5*(gdata.om[q_Bx](i,j,k) +
							 gdata.om[q_Bx](i-1,j,k));
				  B_theta = 0.5*(gdata.om[q_By](i,j,k) +
								 gdata.om[q_By](i,j-1,k));
				  B_phi = 0.5*(gdata.om[q_Bz](i,j,k) +
							   gdata.om[q_Bz](i,j,k-1));
				// } else {
				//   B_r = 0.;
				//   B_theta = 0.;
				//   B_phi = 0.;
				}
//#else
//				REAL B_r(0.), B_theta(0.), B_phi(0.);
//#endif

#if (USE_COROTATION == CRONOS_ON)
				REAL rad = gdata.getCen_x(i);
				REAL v_phiCorot = v_phi - rad*omegaZ*sin(theta);

				// Source for v_r
				nom[q_sx](i,j,k) -=  (rho*(v_phiCorot*v_phi + sqr(v_theta)) +
				                         2.*pres + sqr(B_r))*r_inv;

				// Source for v_theta
				nom[q_sy](i,j,k) -= (-rho*v_r*v_theta + B_r*B_theta +
				                        cot*(rho*v_phiCorot*v_phi + pres +
				                             0.5*(sqr(B_r) + sqr(B_theta) - sqr(B_phi))))*r_inv;
				// Source for v_phi
				nom[q_sz](i,j,k) -= (-rho*v_r*v_phiCorot + B_r*B_phi +
				                        cot*(-rho*v_theta*v_phiCorot + B_phi*B_theta))*r_inv;
#else
				// Source for v_r
				nom[q_sx](i,j,k) -=  (rho*(sqr(v_phi) + sqr(v_theta)) +
				                         2.*pres + sqr(B_r))*r_inv;

				// Source for v_theta
				nom[q_sy](i,j,k) -= (-rho*v_r*v_theta + B_r*B_theta +
				                        cot*(rho*sqr(v_phi) + pres +
				                             0.5*(sqr(B_r) + sqr(B_theta) - sqr(B_phi))))*r_inv;
				// Source for v_phi
				nom[q_sz](i,j,k) -= (-rho*v_r*v_phi + B_r*B_phi +
				                        cot*(-rho*v_phi*v_theta + B_phi*B_theta))*r_inv;
#endif

	
			}
		}
	}

}

#ifdef PHYSDISS
void SourceTerms::src_Geom_Viscosity(Data &gdata, ProblemType &pr,
                                     NumMatrix<REAL,3> nom[N_OMINT], const CronosFluid &fluid) {

	// Set indices:
	int q_rho = fluid.get_q_rho();
	int q_sx = fluid.get_q_sx();
	int q_sy = fluid.get_q_sy();
	int q_sz = fluid.get_q_sz();

	if(gdata.om[q_sx].getName() != "v_x" ||
	   gdata.om[q_sy].getName() != "v_y" ||
	   gdata.om[q_sz].getName() != "v_z") {
		throw CException(" Need characteristic variables for geom sources ");
	}

	Pot & rho   = gdata.om[q_rho];
	Pot & v_rad = gdata.om[q_sx];
	Pot & v_tht = gdata.om[q_sy];
	Pot & v_phi = gdata.om[q_sz];

	for (int k = 0; k <= gdata.mx[2]; ++k){
		for (int j = 0; j <= gdata.mx[1]; ++j){
			REAL theta = gdata.getCen_y(j);
			REAL sin_theta = sin(theta);
			REAL sin_theta_inv = 1./sin_theta;
			REAL cot_theta = 1./tan(theta);
			//    REAL cot_theta_inv = 1./cot_theta;
			for (int i = 0; i <= gdata.mx[0]; ++i){

				REAL r_sph     = gdata.getCen_x(i);
				REAL r_sph_inv = 1./r_sph;

				REAL mu = pr.nu(gdata,i,j,k)*rho(i,j,k);

				REAL h0_inv = 1./gdata.h0(i,j,k);
				REAL h1_inv = 1./gdata.h1(i,j,k);
				REAL h2_inv = 1./gdata.h2(i,j,k);
	

				REAL div_v = ((v_rad(i+1,j,k) - v_rad(i-1,j,k))*gdata.hx[0] +
				              2.*v_rad(i,j,k)*r_sph_inv +
				              (v_tht(i,j+1,k) - v_tht(i,j-1,k))*gdata.hx[1]*h1_inv +
				              (v_phi(i,j,k+1) - v_phi(i,j,k-1))*gdata.hx[2]*h2_inv +
				              v_tht(i,j,k)*r_sph_inv*cot_theta);

				REAL tau_pp = (2.*((v_phi(i,j,k+1) - v_phi(i,j,k-1))*gdata.hx[2]*h2_inv+
				                   v_rad(i,j,k)*r_sph_inv +
				                   v_tht(i,j,k)*r_sph_inv*sin_theta_inv) -
				               twothirds*div_v);

				REAL tau_tt = (2.*((v_tht(i,j+1,k) - 
				                    v_tht(i,j-1,k))*gdata.hx[1]*h1_inv +
				                   v_rad(i,j,k)*r_sph_inv) -
				               twothirds*div_v);

				REAL tau_pr = ((v_rad(i,j,k+1) - v_rad(i,j,k-1))*gdata.hx[2]*h2_inv +
				               (v_phi(i+1,j,k) - v_phi(i-1,j,k))*gdata.hx[0] -
				               v_phi(i,j,k)*r_sph_inv);

				REAL tau_pt = ((v_phi(i,j+1,k) - v_phi(i,j-1,k))*gdata.hx[1]*h1_inv +
				               - v_phi(i,j,k)*r_sph_inv*cot_theta +
				               (v_tht(i,j,k+1) - v_tht(i,j,k-1))*gdata.hx[2]*h2_inv);
		       
				REAL tau_tr = ((v_rad(i,j+1,k) - v_rad(i,j-1,k))*gdata.hx[1]*h1_inv +
				               (v_tht(i+1,j,k) - v_tht(i-1,j,k))*gdata.hx[0]*h0_inv -
				               v_tht(i,j,k)*r_sph_inv);

				nom[q_sx](i,j,k) -= mu*(-tau_pp - tau_tt)*r_sph_inv;
	
				nom[q_sy](i,j,k) -= mu*( tau_pr + tau_pt*cot_theta)*r_sph_inv;

				nom[q_sz](i,j,k) -= mu*( tau_tr - tau_pp*cot_theta)*r_sph_inv;

			}
		}
	}
  
}



void SourceTerms::mod_Geom_Visc_WE(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[N2OMINT], REAL tau[N_OMINT],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	cerr << " sources::mod_Geom_Visc_WE - not adapted for multifluid yet " << endl;
	exit(3);

	REAL r_sph_inv = 1./(gdata.getCen_x(i));
	REAL cot_theta = 1./tan(gdata.getCen_y(j));

	// tau_rr
	tau[1] += -(omWE[Problem.q_sx+dir]*2. + 
	            omWE[Problem.q_sz+dir]*cot_theta)*r_sph_inv*twothirds;

	// tau_rtheta
	tau[2] += -omWE[Problem.q_sy+dir]*r_sph_inv;

	// tau_rphi
	tau[3] += -omWE[Problem.q_sz+dir]*r_sph_inv;

	return;
}


void SourceTerms::mod_Geom_Visc_SN(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[N2OMINT], REAL tau[N_OMINT],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	cerr << " sources::mod_Geom_Visc_SN - not adapted for multifluid yet " << endl;
	exit(3);

	REAL r_sph_inv = 1./(gdata.getCen_x(i));
	REAL cot_theta = 1./tan(gdata.getCen_y(j));

	// tau_phir
	tau[1] += -omWE[Problem.q_sz+dir]*r_sph_inv;

	// tau_phitheta
	tau[2] += -omWE[Problem.q_sz+dir]*r_sph_inv*cot_theta;

	// tau_phiphi
	tau[3] +=  (omWE[Problem.q_sx+dir] + omWE[Problem.q_sy+dir]*2.*cot_theta)*r_sph_inv*twothirds;

	return;
}


void SourceTerms::mod_Geom_Visc_BT(Data &gdata, const ProblemType &Problem,
                                   REAL omWE[N2OMINT], REAL tau[N_OMINT],
                                   const REAL &i, const REAL &j, const REAL &k,
                                   const int &dir) {
	cerr << " sources::mod_Geom_Visc_BT - not adapted for multifluid yet " << endl;
	exit(3);

	REAL r_sph_inv = 1./(gdata.getCen_x(i));
	REAL cot_theta = 1./tan(gdata.getCen_y(j));

	// tau_thetar
	tau[1] += -omWE[Problem.q_sy+dir]*r_sph_inv;

	// tau_thetatheta
	tau[2] +=  (omWE[Problem.q_sx+dir] - 
	            omWE[Problem.q_sy+dir]*cot_theta)*r_sph_inv*twothirds;

	// tau_thetaphi
	tau[3] += -omWE[Problem.q_sz+dir]*r_sph_inv*cot_theta;

	return;
}
#endif
#endif


void SourceTerms::src_Geom(Data &gdata, ProblemType &Problem,
                           NumMatrix<REAL,3> nom[N_OMINT]) {
	//! Appy geometrical source terms for non-Cartesian grids

	// In case of multifluid simulation: loop over all fluids:
#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	int numFluids = gdata.fluids->get_numFluids();
	for (int iFluid=0; iFluid<numFluids; ++iFluid) {
		CronosFluid fluid = gdata.fluids->fluids[iFluid];
#else
		CronosFluid fluid = gdata.fluid;
#endif

	src_Geom_Hydro(gdata, Problem, nom, fluid);
#ifdef PHYSDISS
	src_Geom_Viscosity(gdata, Problem, nom, fluid);
#endif

#if (GEOM == CYLINDRICAL)
	//	if(gdata.axis_treatment(gdata.phi_dir) == 1) {
	if(gdata.get_singularity_treatment(0) == 1) {
		src_Axis(gdata, nom, fluid);
	}
#endif

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	}
#endif

}


void SourceTerms::src_Axis(Data &gdata, NumMatrix<REAL,3> nom[N_OMINT], const CronosFluid &fluid) {

	int q_Bx = fluid.get_q_Bx();
	int q_By = fluid.get_q_By();
	int q_Bz = fluid.get_q_Bz();
	int fluid_type = fluid.get_fluid_type();

#if (FLUID_TYPE == CRONOS_MHD)
#if (GEOM == 2)
	if(IntegrateA) { 
		for (int k = 0; k <= gdata.mx[2]; k++) {
			for (int j = 0; j <= gdata.mx[1]; j++) {
				nom[q_By](0,j,k) = 0.;
				// nom[q_Bz](0,j,k) = 0.;
			}
		}
	} else {
		for (int k = 0; k <= gdata.mx[2]; k++) {
			for (int j = 0; j <= gdata.mx[1]; j++) {
				nom[q_Bx](-1,j,k) = 0.;
			}
		}
	}
#endif
#endif

}


void SourceTerms::src_rhs(Data &gdata, ProblemType &Problem,
                          NumMatrix<REAL,3> nom[N_OMINT]) {
#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	int numFluids = gdata.fluids->get_numFluids();
	for (int iFluid=0; iFluid<numFluids; ++iFluid) {
		CronosFluid fluid = gdata.fluids->fluids[iFluid];
#else
		CronosFluid fluid = gdata.fluid;
#endif

	src_rhs_Ideal(gdata, Problem, nom, fluid);
#ifdef PHYSDISS
	src_rhs_Viscosity(gdata, Problem, nom, fluid);
#endif

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	}
#endif
}


void SourceTerms::src_rhs_Ideal(Data &gdata, ProblemType &pr,
                                NumMatrix<REAL,3> nom[N_OMINT], const CronosFluid &fluid) {

	// Source terms as, e.g., for RHS of thermal energy equations

	// Set indices:
	int q_sx = fluid.get_q_sx();
	int q_sy = fluid.get_q_sy();
	int q_sz = fluid.get_q_sz();
	int q_Eadd = fluid.get_q_Eadd();

	int fluid_type = fluid.get_fluid_type();

	if(gdata.om[q_sx].getName() != "v_x" || gdata.om[q_sy].getName() != "v_y" ||
	   gdata.om[q_sz].getName() != "v_z") {
		throw CException(" Need characteristic variables for geom sources ");
	}

	if(gdata.om[q_Eadd].getName() == "Etherm") {
		for (int k = 0; k <= gdata.mx[2]; ++k) {
			for (int j = 0; j <= gdata.mx[1]; ++j) {
				for (int i = 0; i <= gdata.mx[0]; ++i) {
	  
					REAL divV = ((gdata.om[q_sx](i+1,j,k) -
					              gdata.om[q_sx](i-1,j,k))*gdata.hx[0] +
					             (gdata.om[q_sy](i,j+1,k) -
					              gdata.om[q_sy](i,j-1,k))*gdata.hx[1] +
					             (gdata.om[q_sz](i,j,k+1) -
					              gdata.om[q_sz](i,j,k-1))*gdata.hx[2]);
	  
					nom[q_Eadd](i,j,k) += (pr.gamma - 1.)*gdata.om[q_Eadd](i,j,k)*divV;

				}
			}
		}
	}

}


#ifdef PHYSDISS
void SourceTerms::src_rhs_Viscosity(Data &gdata, ProblemType &pr,
                                    NumMatrix<REAL,3> nom[N_OMINT], const CronosFluid &fluid) {

	// Set indices:
	int q_rho = fluid.get_q_rho();
	int q_sx = fluid.get_q_sx();
	int q_sy = fluid.get_q_sy();
	int q_sz = fluid.get_q_sz();
	int q_Bx = fluid.get_q_Bx();
	int q_By = fluid.get_q_By();
	int q_Bz = fluid.get_q_Bz();
	int q_Eadd = fluid.get_q_Eadd();


	if(gdata.om[q_sx].getName() != "v_x" ||
	   gdata.om[q_sy].getName() != "v_y" ||
	   gdata.om[q_sz].getName() != "v_z") {
		throw CException(" Need characteristic variables for geom sources ");
	}


	// Viscosity source terms for thermal energy / entropy
	for (int k = 0; k <= gdata.mx[2]; ++k) {
		for (int j = 0; j <= gdata.mx[1]; ++j) {
			for (int i = 0; i <= gdata.mx[0]; ++i) {

				REAL dxux = (gdata.om[q_sx](i+1,j,k) -
				             gdata.om[q_sx](i-1,j,k))*gdata.hx[0];
				REAL dyux = (gdata.om[q_sx](i,j+1,k) -
				             gdata.om[q_sx](i,j-1,k))*gdata.hx[1];
				REAL dzux = (gdata.om[q_sx](i,j,k+1) -
				             gdata.om[q_sx](i,j,k-1))*gdata.hx[2];
				REAL dxuy = (gdata.om[q_sy](i+1,j,k) -
				             gdata.om[q_sy](i-1,j,k))*gdata.hx[0];
				REAL dyuy = (gdata.om[q_sy](i,j+1,k) -
				             gdata.om[q_sy](i,j-1,k))*gdata.hx[1];
				REAL dzuy = (gdata.om[q_sy](i,j,k+1) -
				             gdata.om[q_sy](i,j,k-1))*gdata.hx[2];
				REAL dxuz = (gdata.om[q_sz](i+1,j,k) -
				             gdata.om[q_sz](i-1,j,k))*gdata.hx[0];
				REAL dyuz = (gdata.om[q_sz](i,j+1,k) -
				             gdata.om[q_sz](i,j-1,k))*gdata.hx[1];
				REAL dzuz = (gdata.om[q_sz](i,j,k+1) -
				             gdata.om[q_sz](i,j,k-1))*gdata.hx[2];
				REAL divv = dxux + dyuy + dzuz;

				REAL tau11 = 2.*dxux - divv*twothirds;
				REAL tau12 = dxuy + dyux;
				REAL tau13 = dxuz + dzux;
				REAL tau22 = 2.*dyuy - divv*twothirds;
				REAL tau23 = dyuz + dzuy;
				REAL tau33 = 2.*dzuz - divv*twothirds;
        
				REAL rho = gdata.om[q_rho](i,j,k);

#ifdef ENTROPY
				REAL fac = (pr.gamma - 1.)*pow(rho,2.-pr.gamma)*pr.nu(gdata,i,j,k);
#else
				REAL fac = pr.nu(gdata,i,j,k)*rho;
#endif
	
				nom[q_Eadd](i,j,k) -= fac*(tau11*dxux + tau12*dyux +
				                              tau13*dzux + tau12*dxuy + 
				                              tau22*dyuy + tau23*dzuy +
				                              tau13*dxuz + tau23*dyuz + 
				                              tau33*dzuz);
			}
		}
	}  


	// Resistivity source terms for thermal energy / entropy

	for (int k = 0; k <= gdata.mx[2]; ++k) {
		for (int j = 0; j <= gdata.mx[1]; ++j) {
			for (int i = 0; i <= gdata.mx[0]; ++i) {

				REAL eta = pr.eta(gdata, i,j,k);
        
				double dyBxM = (gdata.om[q_Bx](i-1,j+1,k) -
				                gdata.om[q_Bx](i-1,j-1,k))*gdata.hx[1];
				double dyBxP = (gdata.om[q_Bx](i  ,j+1,k) -
				                gdata.om[q_Bx](i  ,j-1,k))*gdata.hx[1];
				double dzBxM = (gdata.om[q_Bx](i-1,j,k+1) -
				                gdata.om[q_Bx](i-1,j,k-1))*gdata.hx[2];
				double dzBxP = (gdata.om[q_Bx](i  ,j,k+1) -
				                gdata.om[q_Bx](i  ,j,k-1))*gdata.hx[2];
        
				double dxByM = (gdata.om[q_By](i+1,j  ,k) -
				                gdata.om[q_By](i-1,j  ,k))*gdata.hx[0];
				double dxByP = (gdata.om[q_By](i+1,j-1,k) -
				                gdata.om[q_By](i-1,j-1,k))*gdata.hx[0];
				double dzByM = (gdata.om[q_By](i,j  ,k+1) -
				                gdata.om[q_By](i,j  ,k-1))*gdata.hx[2];
				double dzByP = (gdata.om[q_By](i,j-1,k+1) -
				                gdata.om[q_By](i,j-1,k-1))*gdata.hx[2];
        
				double dxBzP = (gdata.om[q_Bz](i+1,j,k  ) -
				                gdata.om[q_Bz](i-1,j,k  ))*gdata.hx[0];
				double dxBzM = (gdata.om[q_Bz](i+1,j,k-1) -
				                gdata.om[q_Bz](i-1,j,k-1))*gdata.hx[0];
				double dyBzP = (gdata.om[q_Bz](i,j+1,k  ) -
				                gdata.om[q_Bz](i,j-1,k  ))*gdata.hx[1];
				double dyBzM = (gdata.om[q_Bz](i,j+1,k-1) -
				                gdata.om[q_Bz](i,j-1,k-1))*gdata.hx[1];
        
				double JxP = dyBzP - dzByP;
				double JxM = dyBzM - dzByM;
				double JyP = dzBxP - dxByP;
				double JyM = dzBxM - dxByM;
				double JzP = dxByP - dyBxP;
				double JzM = dxByM - dyBxM;
				double Jxsqr = (sqr(JxP) + sqr(JxM) + JxP*JxM)*onethird;
				double Jysqr = (sqr(JyP) + sqr(JyM) + JyP*JyM)*onethird;
				double Jzsqr = (sqr(JzP) + sqr(JzM) + JzP*JzM)*onethird;

#ifdef ENTROPY
				nom[q_Eadd](i,j,k) -= eta*(gamma-1.)(Jxsqr + Jysqr + Jzsqr)*pow(gdata.om[q_rho](i,j,k),1-gamma);
#else
				nom[q_Eadd](i,j,k) -= eta*(Jxsqr + Jysqr + Jzsqr);
#endif
			}
		}
	}

}
#endif
