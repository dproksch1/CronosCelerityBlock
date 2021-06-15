#include "gridgen.H"
#include "vector.H"
#include "CException.H"
#include <string.h>
#include <iomanip>
#include <cmath>

REAL HyperbolicSolver::singlestep(Data &gdata, gridFunc &gfunc,
                                  ProblemType &Problem, int n)
{

	// Typical definition of variables for MHD run:
	//   om[0] : rho                 rho
	//   om[1] : rho v_x             rvx
	//   om[2] : rho v_y             rvy
	//   om[3] : rho v_z             rvz
	//   om[4] : B_x                 bx
	//   om[5] : B_y                 by
	//   om[6] : B_z                 bz
	//   om[7] : e_th / e_ges        e (if N_OMINT == 8)
  

#ifdef PHYSDISS
	double comp_max_glob = 0.;
#endif

	const int N2OMINT = 2*N_OMINT;
	const int N_VAR = N_OMINT-4;
	const double norm = 1./sqrt(2.); // Some constant for velocities


	int q, i, j, k, kold;        // q numbers vars, {i,j,k} for grid pts


	// Check if there are negative Temperatures
	if(gdata.rank == 0) {
		if(N_OMINT >= 8 && TempNeg !=0) {
			cout << " Negative Temperatures: " << TempNeg << endl;
			TempNeg = 0;
		}
	}

	/* Definition of temporary fields holding the change of the
	   individual variables - field is deleted when not needed anymore
	   mx[]+1 = grid extension in [x,y,z]
	*/
  
	int lbound[3],ubound[3];

	for(int i=0;i<3;++i){
		lbound[i]=0;
		ubound[i]=gdata.mx[i];
	}
	NumMatrix<REAL,3>* nom;

	nom = new NumMatrix<REAL,3> [N_OMINT];
	for (q = 0; q < N_OMINT; ++q) {
		nom[q].resize(lbound,ubound); 
		nom[q].clear();
	}
	if(IntegrateA) {
		nom[4].resize(Index::set(0,0,0),
		              Index::set(gdata.mx[0],gdata.mx[1]+1,gdata.mx[2]+1));
		nom[5].resize(Index::set(0,0,0),
		              Index::set(gdata.mx[0]+1,gdata.mx[1],gdata.mx[2]+1));
		nom[6].resize(Index::set(0,0,0),
		              Index::set(gdata.mx[0]+1,gdata.mx[1]+1,gdata.mx[2]));
	} else {
		nom[4].resize(Index::set(-1,0,0),
		              Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
		nom[5].resize(Index::set(0,-1,0),
		              Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
		nom[6].resize(Index::set(0,0,-1),
		              Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
	}


	nom[4].clear();
	nom[5].clear();
	nom[6].clear();

	for(int i=0;i<3;++i){
		lbound[i]=-1;
		ubound[i]=gdata.mx[i]+1;
	}
	NumMatrix<REAL,3> HXpmWE[N_OMINT], HXpmSN[N_OMINT], HXpmBT[N_OMINT];
	for (q = 0; q < N_OMINT; ++q) {
		if(q>=4 && q<=6) {
			HXpmWE[q].resize(Index::set(-1,-1,-1),Index::set(-1,-1,-1));
			HXpmSN[q].resize(Index::set(-1,-1,-1),Index::set(-1,-1,-1));
			HXpmBT[q].resize(Index::set(-1,-1,-1),Index::set(-1,-1,-1));
		} else {
			HXpmWE[q].resize(lbound,ubound);
			HXpmSN[q].resize(lbound,ubound);
			HXpmBT[q].resize(lbound,ubound);
		}
		HXpmWE[q].clear();
		HXpmSN[q].clear();
		HXpmBT[q].clear();
	}


#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	if(N_OMINT > 7) {
		Trafo->computeEntropyFromE(gdata, gfunc, Problem);
	}
#endif

	if(N_OMINT > 7) {
		if(gdata.om[7].getName() == "Etherm") {
			Trafo->TransEth2E(gdata, gfunc, Problem);
		} else if(gdata.om[7].getName() == "Temp") {
			Trafo->TransT2E(gdata, gfunc, Problem);
		}
	}

#ifdef ANGULAR_MOMENTUM
	Trafo->TransVel2AngMom(gdata, gfunc, Problem);
#else
	Trafo->TransVel2Momen(gdata, gfunc, Problem);
#endif

	// Saving old variables for Runge-Kutta steps:
	if (n == 0) {
		for (q = 0; q < N_OMINT; ++q){
			if((q < 4 || q > 6) || !IntegrateA) {
				gdata.om[N_OMEGA+q] = gdata.om[q];         // save om_n
			} else {
				gdata.om[N_OMEGA+q] = gdata.om[N_OMINT+N_ADD+q-4];
			}
		}
	}

	if(eos == NULL) {
		eos = new EquationOfState(Problem);
	}
	if(sources == NULL) {
		sources = new SourceTerms(gdata, Problem, IntegrateA);
	}


	REAL pres;
	REAL v_max_m[3], v_cor_m[3];
	REAL v_max_p[3], v_cor_p[3];
	REAL ommx[N_OMINT+3],flmx[N_OMINT];
	//   NumMatrix<REAL,1> dissWE[N_OMINT], dissSN[N_OMINT], dissBT[N_OMINT];
	//   NumMatrix<REAL,1> compWE[N_OMINT] ,compSN[N_OMINT], compBT[N_OMINT];
	//   NumMatrix<REAL,1> dissmy[N_OMINT], compmy[N_OMINT];
	//   NumMatrix<REAL,2> dissmz[N_OMINT], compmz[N_OMINT];
	NumMatrix<REAL,1> ommy[N_OMINT+3],flmy[N_OMINT+1];
	NumMatrix<REAL,1> vzm_old, vzp_old;
	NumMatrix<REAL,2> ommz[N_OMINT+3],flmz[N_OMINT+3];;
	NumMatrix<REAL,2> vym_old, vyp_old, vxm_old, vxp_old;
	NumMatrix<REAL,3> Ex, Ey, Ez;
	NumMatrix<REAL,3> A, Ap, Am;
  

	vxm_old.resize(Index::set(-2,-2),Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	vxp_old.resize(Index::set(-2,-2),Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	vym_old.resize(Index::set(-2,-2),Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	vyp_old.resize(Index::set(-2,-2),Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	vzm_old.resize(Index::set(-2),Index::set(gdata.mx[0]+2));
	vzp_old.resize(Index::set(-2),Index::set(gdata.mx[0]+2));

	for (q = 0; q < N_OMINT+3; ++q) {
		ommy[q].resize(Index::set(-2),Index::set(gdata.mx[0]+2));
		ommz[q].resize(Index::set(-2,-2),Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	}

	for (q = 0; q < N_OMINT; ++q) {
		flmy[q].resize(Index::set(-2),Index::set(gdata.mx[0]+2));
		flmy[q].clear();
		flmz[q].resize(Index::set(-2,-2),Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	}

	flmy[N_OMINT].resize(Index::set(-2),Index::set(gdata.mx[0]+2));

	flmz[N_OMINT].resize(Index::set(-2,-2),
	                     Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	flmz[N_OMINT+1].resize(Index::set(-2,-2),
	                       Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	flmz[N_OMINT+2].resize(Index::set(-2,-2),
	                       Index::set(gdata.mx[0]+2,gdata.mx[1]+2));
	A.resize(Index::set(0,0,0),Index::set( 1, 1, 1));
	Am.resize(Index::set(0,0,0),Index::set( 1, 1, 1));
	Ap.resize(Index::set(0,0,0),Index::set( 1, 1, 1));

	Ex.resize(Index::set(-1,-1,-1),
	          Index::set(gdata.mx[0]+1,gdata.mx[1]+1,gdata.mx[2]+1));
	Ey.resize(Index::set(-1,-1,-1),
	          Index::set(gdata.mx[0]+1,gdata.mx[1]+1,gdata.mx[2]+1));
	Ez.resize(Index::set(-1,-1,-1),
	          Index::set(gdata.mx[0]+1,gdata.mx[1]+1,gdata.mx[2]+1));

  
	REAL omWE[N2OMINT], omSN[N2OMINT], omBT[N2OMINT];
	REAL omWS[N_OMINT], omES[N_OMINT], omWN[N_OMINT], omEN[N_OMINT];
	REAL omWB[N_OMINT], omEB[N_OMINT], omWT[N_OMINT], omET[N_OMINT];
	REAL omSB[N_OMINT], omNB[N_OMINT], omST[N_OMINT], omNT[N_OMINT];
	REAL ExSB(0.), ExNB(0.), ExST(0.), ExNT(0.), EyWB(0.), EyEB(0.);
	REAL EyWT(0.), EyET(0.), EzWS(0.), EzES(0.), EzWN(0.), EzEN(0.);

	REAL fluxWE[N2OMINT], fluxSN[N2OMINT], fluxBT[N2OMINT];
	REAL domdx[N2OMINT], domdy[N2OMINT], domdz[N2OMINT];

	for(int i=0; i<N2OMINT; ++i){
		fluxSN[i] = 0.;
		omWE[i] = 0.;
		omSN[i] = 0.;
		omBT[i] = 0.;
	}

	REAL cfl_lin(0.);
#ifdef PHYSDISS
	REAL cfl_eta(0.), cfl_eta_loc(0.);
#endif

	// ---------------------------------------------------------------	      
	//      Using primitive variables for reconstruction
	//----------------------------------------------------------------

	// As soon as the variables for the previous timestep are saved, we
	// transform from vector potential to magnetic induction

	//  TransPot2Mag(gdata, gfunc, Problem);

#ifdef ANGULAR_MOMENTUM
	Trafo->TransAngMom2Vel(gdata, gfunc, Problem);
#else
	Trafo->TransMomen2Vel(gdata, gfunc, Problem);
#endif

	if(N_OMINT >= 8) {
		if(thermal) {
			Trafo->TransE2Eth(gdata, gfunc, Problem);
		} else {
			Trafo->TransE2T(gdata, gfunc, Problem);
		}
	}

	gettimeofday(&tick, 0);
	cstart = clock();
  
	gdata.CheckNan(1);

	double dxinv[3];
	dxinv[0] = 1./gdata.dx[0];
	dxinv[1] = 1./gdata.dx[1];
	dxinv[2] = 1./gdata.dx[2];


	kold = -2;
	int jold(-2);
	for (k = -2; k <= gdata.mx[2]+1; ++k){
		if(k > kold) {
			kold = k;
			for (j = -2; j <= gdata.mx[1]+1; ++j){
				for (i = -2; i <= gdata.mx[0]+1; ++i){
					flmz[N_OMINT+1](i,j) = flmz[N_OMINT](i,j);
					flmz[N_OMINT+2](i,j) = flmz[6](i,j);
				}
			}
		}
		for (j = -2; j <= gdata.mx[1]+1; ++j){
			if(j > jold) {
				for (i = -2; i <= gdata.mx[0]+1; ++i){
					flmy[N_OMINT](i)     = flmy[5](i);
				}
			}
			for (i = -2; i <= gdata.mx[0]+1; ++i){

				// Compute shift of collocation point:
				REAL shiftWE(0.), shiftSN(0.), shiftBT(0.);
#ifdef GEOM
#if GEOM != CARTESIAN
				shiftWE = sources->shift_Geom_WE(gdata, i,j,k);
				shiftSN = sources->shift_Geom_SN(gdata, i,j,k);
				shiftBT = sources->shift_Geom_BT(gdata, i,j,k);
#endif
#endif

				for (q = 0; q < 4; ++q){
					// Reconstruction for density and velocity

					// left-handed, central and right-handed gradients for all
					// directions
					REAL omxp =  gdata.om[q](i+1,j  ,k  ) - gdata.om[q](i  ,j  ,k  );
					REAL omx0 = (gdata.om[q](i+1,j  ,k  ) - gdata.om[q](i-1,j  ,k  ))/2.;
					REAL omxm =  gdata.om[q](i  ,j  ,k  ) - gdata.om[q](i-1,j  ,k  );
					REAL omyp =  gdata.om[q](i  ,j+1,k  ) - gdata.om[q](i  ,j  ,k  );
					REAL omy0 = (gdata.om[q](i  ,j+1,k  ) - gdata.om[q](i  ,j-1,k  ))/2.;
					REAL omym =  gdata.om[q](i  ,j  ,k  ) - gdata.om[q](i  ,j-1,k  );
					REAL omzp =  gdata.om[q](i  ,j  ,k+1) - gdata.om[q](i  ,j  ,k  );
					REAL omz0 = (gdata.om[q](i  ,j  ,k+1) - gdata.om[q](i  ,j  ,k-1))/2.;
					REAL omzm =  gdata.om[q](i  ,j  ,k  ) - gdata.om[q](i  ,j  ,k-1);
	  
	  
					A(0,0,0) =  gdata.om[q](i,j,k);
	  
					// gradient in x-direction
					//	    A(1,0,0) = Limiter.compute(omxp,omx0,omxm);
					A(1,0,0) = Limiter->compute(omxp,omx0,omxm);
					// gradient in y-direction
					A(0,1,0) = Limiter->compute(omyp,omy0,omym);
					// gradient in z-direction
					A(0,0,1) = Limiter->compute(omzp,omz0,omzm);
	  
					// Here N_OMINT signifies the values at the upper boundary
					// point values on lower (higher) x-boundary
					omWE[q] = A(0,0,0) - (0.5+shiftWE)*A(1,0,0);
	  
					omWE[q+N_OMINT] = omWE[q] + A(1,0,0);

					domdx[q] = A(1,0,0);
					// point values on lower (higher) y-boundary
					omSN[q] = A(0,0,0) - (0.5+shiftSN)*A(0,1,0);

					omSN[q+N_OMINT] = omSN[q] + A(0,1,0);
					domdy[q] = A(0,1,0);
					// point values on lower (higher) z-boundary
					omBT[q] = A(0,0,0) - (0.5+shiftBT)*A(0,0,1);
					omBT[q+N_OMINT] = omBT[q] + A(0,0,1);
					domdz[q] = A(0,0,1);

	  
					// Reconstruction of corner values for computation of emf
					if(Problem.mag) {
						if(q > 0) {
							omWS[q] = A(0,0,0) - ((0.5+shiftWE)*A(1,0,0) +
							                      (0.5+shiftSN)*A(0,1,0));
							omES[q] = omWS[q] + A(1,0,0);
							omWN[q] = omWS[q] + A(0,1,0);
							omEN[q] = omES[q] + A(0,1,0);
	      
							omWB[q] = A(0,0,0) - ((0.5+shiftWE)*A(1,0,0) +
							                      (0.5+shiftBT)*A(0,0,1));
							omEB[q] = omWB[q] + A(1,0,0);
							omWT[q] = omWB[q] + A(0,0,1);
							omET[q] = omEB[q] + A(0,0,1);
	      
							omSB[q] = A(0,0,0) - ((0.5+shiftSN)*A(0,1,0) +
							                      (0.5+shiftBT)*A(0,0,1));
							omNB[q] = omSB[q] + A(0,1,0);
							omST[q] = omSB[q] + A(0,0,1);
							omNT[q] = omNB[q] + A(0,0,1);

						}
					}
				}
	

				if(Problem.mag){
					/* 	Reconstruction for magnetic field: individually for
						all directions and components of magnetic field.  For
						this keep in mind that magnetic field is defined on
						the following locations
						Bx(i,j,k) = om[4](i,j,k) <-> Bx(i+1/2,j,k)
						By(i,j,k) = om[5](i,j,k) <-> By(i,j+1/2,k)
						Bz(i,j,k) = om[6](i,j,k) <-> Bz(i,j,k+1/2)
					*/

					// Computing all relevant local derivatives for Bx
					REAL dyBxpr =  gdata.om[4](i  ,j+1,k  ) - gdata.om[4](i  ,j  ,k  );
					REAL dyBx0r = (gdata.om[4](i  ,j+1,k  ) -
					               gdata.om[4](i  ,j-1,k  ))/2.;
					REAL dyBxmr =  gdata.om[4](i  ,j  ,k  ) - gdata.om[4](i  ,j-1,k  );
					REAL dzBxpr =  gdata.om[4](i  ,j  ,k+1) - gdata.om[4](i  ,j  ,k  );
					REAL dzBx0r = (gdata.om[4](i  ,j  ,k+1) -
					               gdata.om[4](i  ,j  ,k-1))/2.;
					REAL dzBxmr =  gdata.om[4](i  ,j  ,k  ) - gdata.om[4](i  ,j  ,k-1);
					REAL dyBxpl =  gdata.om[4](i-1,j+1,k  ) - gdata.om[4](i-1,j  ,k  );
					REAL dyBx0l = (gdata.om[4](i-1,j+1,k  ) -
					               gdata.om[4](i-1,j-1,k  ))/2.;
					REAL dyBxml =  gdata.om[4](i-1,j  ,k  ) - gdata.om[4](i-1,j-1,k  );
					REAL dzBxpl =  gdata.om[4](i-1,j  ,k+1) - gdata.om[4](i-1,j  ,k  );
					REAL dzBx0l = (gdata.om[4](i-1,j  ,k+1) -
					               gdata.om[4](i-1,j  ,k-1))/2.;
					REAL dzBxml =  gdata.om[4](i-1,j  ,k  ) - gdata.om[4](i-1,j  ,k-1);
	  
					// Computing all relevant local derivatives for By
					REAL dxBypr =  gdata.om[5](i+1,j  ,k  ) - gdata.om[5](i  ,j  ,k  );
					REAL dxBy0r = (gdata.om[5](i+1,j  ,k  ) -
					               gdata.om[5](i-1,j  ,k  ))/2.;
					REAL dxBymr =  gdata.om[5](i  ,j  ,k  ) - gdata.om[5](i-1,j  ,k  );
					REAL dzBypr =  gdata.om[5](i  ,j  ,k+1) - gdata.om[5](i  ,j  ,k  );
					REAL dzBy0r = (gdata.om[5](i  ,j  ,k+1) -
					               gdata.om[5](i  ,j  ,k-1))/2.;
					REAL dzBymr =  gdata.om[5](i  ,j  ,k  ) - gdata.om[5](i  ,j  ,k-1);
					REAL dxBypl =  gdata.om[5](i+1,j-1,k  ) - gdata.om[5](i  ,j-1,k  );
					REAL dxBy0l = (gdata.om[5](i+1,j-1,k  ) -
					               gdata.om[5](i-1,j-1,k  ))/2.;
					REAL dxByml =  gdata.om[5](i  ,j-1,k  ) - gdata.om[5](i-1,j-1,k  );
					REAL dzBypl =  gdata.om[5](i  ,j-1,k+1) - gdata.om[5](i  ,j-1,k  );
					REAL dzBy0l = (gdata.om[5](i  ,j-1,k+1) - 
					               gdata.om[5](i  ,j-1,k-1))/2.;
					REAL dzByml =  gdata.om[5](i  ,j-1,k  ) - gdata.om[5](i  ,j-1,k-1);
	  
					// Computing all relevant local derivatives for Bz
					REAL dxBzpr =  gdata.om[6](i+1,j  ,k  ) - gdata.om[6](i  ,j  ,k  );
					REAL dxBz0r = (gdata.om[6](i+1,j  ,k  ) - 
					               gdata.om[6](i-1,j  ,k  ))/2.;
					REAL dxBzmr =  gdata.om[6](i  ,j  ,k  ) - gdata.om[6](i-1,j  ,k  );
					REAL dyBzpr =  gdata.om[6](i  ,j+1,k  ) - gdata.om[6](i  ,j  ,k  );
					REAL dyBz0r = (gdata.om[6](i  ,j+1,k  ) - 
					               gdata.om[6](i  ,j-1,k  ))/2.;
					REAL dyBzmr =  gdata.om[6](i  ,j  ,k  ) - gdata.om[6](i  ,j-1,k  );
					REAL dxBzpl =  gdata.om[6](i+1,j  ,k-1) - gdata.om[6](i  ,j  ,k-1);
					REAL dxBz0l = (gdata.om[6](i+1,j  ,k-1) - 
					               gdata.om[6](i-1,j  ,k-1))/2.;
					REAL dxBzml =  gdata.om[6](i  ,j  ,k-1) - gdata.om[6](i-1,j  ,k-1);
					REAL dyBzpl =  gdata.om[6](i  ,j+1,k-1) - gdata.om[6](i  ,j  ,k-1);
					REAL dyBz0l = (gdata.om[6](i  ,j+1,k-1) - 
					               gdata.om[6](i  ,j-1,k-1))/2.;
					REAL dyBzml =  gdata.om[6](i  ,j  ,k-1) - gdata.om[6](i  ,j-1,k-1);
	  

					Am(0,1,0) = Limiter->compute(dyBxpl,dyBx0l,dyBxml);
					Ap(0,1,0) = Limiter->compute(dyBxpr,dyBx0r,dyBxmr);
					Am(0,0,1) = Limiter->compute(dzBxpl,dzBx0l,dzBxml);
					Ap(0,0,1) = Limiter->compute(dzBxpr,dzBx0r,dzBxmr);

					A(0,1,0) = (Am(0,1,0) + Ap(0,1,0))*0.5;
					A(0,0,1) = (Am(0,0,1) + Ap(0,0,1))*0.5;

					// Saving gradients for dissipation:
					domdy[4]         = Am(0,1,0);
					domdy[4+N_OMINT] = Ap(0,1,0);
					domdz[4]         = Am(0,0,1);
					domdz[4+N_OMINT] = Ap(0,0,1);

					// Computing point values for Bx in x-direction
					omWE[4]         =  gdata.om[4](i-1,j,k);
					omWE[4+N_OMINT] =  gdata.om[4](i  ,j,k);

					// Computing point values for Bx in y-direction
					omSN[4]         = ((gdata.om[4](i-1,j,k) + gdata.om[4](i,j,k))*0.5 -
					                   A(0,1,0)*(0.5+shiftSN));;
					omSN[4+N_OMINT] = (omSN[4] + A(0,1,0));

					// Computing point values for Bx in z-direction
					omBT[4]         = ((gdata.om[4](i-1,j,k) + gdata.om[4](i,j,k))*0.5 -
					                   A(0,0,1)*(0.5+shiftBT));
					omBT[4+N_OMINT] = (omBT[4] + A(0,0,1));

					// Computing point values for Bx at xy-corners
					omWS[4] = gdata.om[4](i-1,j,k) - Am(0,1,0)*(0.5+shiftSN);
					omWN[4] = omWS[4] + Am(0,1,0);
					omES[4] = gdata.om[4](i  ,j,k) - Ap(0,1,0)*(0.5+shiftSN);
					omEN[4] = omES[4] + Ap(0,1,0);


					// Computing point values for Bx at xz-corners
					omWB[4] = gdata.om[4](i-1,j,k) - Am(0,0,1)*(0.5+shiftBT);
					omWT[4] = omWB[4] + Am(0,0,1);
					omEB[4] = gdata.om[4](i  ,j,k) - Ap(0,0,1)*(0.5+shiftBT);
					omET[4] = omEB[4] + Ap(0,0,1);


					Am(1,0,0) = Limiter->compute(dxBypl,dxBy0l,dxByml);
					Ap(1,0,0) = Limiter->compute(dxBypr,dxBy0r,dxBymr);
					Am(0,0,1) = Limiter->compute(dzBypl,dzBy0l,dzByml);
					Ap(0,0,1) = Limiter->compute(dzBypr,dzBy0r,dzBymr);

					A(1,0,0) = (Am(1,0,0) + Ap(1,0,0))*0.5;
					A(0,0,1) = (Am(0,0,1) + Ap(0,0,1))*0.5;
	  
					// Saving gradients for dissipation:
					domdx[5]         = Am(1,0,0);
					domdx[5+N_OMINT] = Ap(1,0,0);
					domdz[5]         = Am(0,0,1);
					domdz[5+N_OMINT] = Ap(0,0,1);

					// Computing point values for By in x-direction
					omWE[5]         = ((gdata.om[5](i,j-1,k) + gdata.om[5](i,j,k))*0.5 -
					                   A(1,0,0)*(0.5+shiftWE));
					omWE[5+N_OMINT] = (omWE[5] + A(1,0,0));

					// Computing point values for By in y-direction
					omSN[5]         =  gdata.om[5](i,j-1,k);
					omSN[5+N_OMINT] =  gdata.om[5](i,j  ,k);
	  
					// Computing point values for By in z-direction
					omBT[5]         = ((gdata.om[5](i,j-1,k) + gdata.om[5](i,j,k))*0.5 -
					                   A(0,0,1)*(0.5+shiftBT));
					omBT[5+N_OMINT] = (omBT[5] + A(0,0,1));

					// Computing point values for By at xy-corners
					omWS[5] = gdata.om[5](i,j-1,k) - 0.5*Am(1,0,0);
					omES[5] = omWS[5] + Am(1,0,0);
					omWN[5] = gdata.om[5](i,j,k) - 0.5*Ap(1,0,0);
					omEN[5] = omWN[5] + Ap(1,0,0);

					// Computing point values for By at yz-corners
					omSB[5] = gdata.om[5](i,j-1,k) - 0.5*Am(0,0,1);
					omST[5] = omSB[5] + Am(0,0,1);
					omNB[5] = gdata.om[5](i,j,k) - 0.5*Ap(0,0,1);
					omNT[5] = omNB[5] + Ap(0,0,1);

	  
					Am(1,0,0) = Limiter->compute(dxBzpl,dxBz0l,dxBzml);
					Ap(1,0,0) = Limiter->compute(dxBzpr,dxBz0r,dxBzmr);
					Am(0,1,0) = Limiter->compute(dyBzpl,dyBz0l,dyBzml);
					Ap(0,1,0) = Limiter->compute(dyBzpr,dyBz0r,dyBzmr);

					A(1,0,0) = (Am(1,0,0) + Ap(1,0,0))*0.5;
					A(0,1,0) = (Am(0,1,0) + Ap(0,1,0))*0.5;

					// Saving gradients for dissipation:
					domdx[5]         = Am(1,0,0);
					domdx[5+N_OMINT] = Ap(1,0,0);
					domdy[5]         = Am(0,1,0);
					domdy[5+N_OMINT] = Ap(0,1,0);

					// Computing point values for Bz in x-direction
					omWE[6]         = (gdata.om[6](i,j,k-1) + 
					                   gdata.om[6](i,j,k  ) - A(1,0,0))*0.5;
					omWE[6+N_OMINT] = (omWE[6] + A(1,0,0));

					// Computing point values for Bz in y-direction
					omSN[6]         = (gdata.om[6](i,j,k-1) + 
					                   gdata.om[6](i,j,k  ) - A(0,1,0))*0.5;
					omSN[6+N_OMINT] = (omSN[6] + A(0,1,0));

					// Computing point values for Bz in z-direction
					omBT[6]         = gdata.om[6](i,j,k-1);
					omBT[6+N_OMINT] = gdata.om[6](i,j,k  );

					// Computing point values for Bz at xz-corners
					omWB[6] = gdata.om[6](i,j,k-1) - 0.5*Am(1,0,0);
					omEB[6] = omWB[6] + Am(1,0,0);
					omWT[6] = gdata.om[6](i,j,k) - 0.5*Ap(1,0,0);
					omET[6] = omWT[6] + Ap(1,0,0);


					// Computing point values for Bz at yz-corners
					omSB[6] = gdata.om[6](i,j,k-1) - 0.5*Am(0,1,0);
					omNB[6] = omSB[6] + Am(0,1,0);
					omST[6] = gdata.om[6](i,j,k) - 0.5*Ap(0,1,0);
					omNT[6] = omST[6] + Ap(0,1,0);


				}


				for (q = 7; q < N_OMINT; ++q){
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
					if(N_OMINT > 8 && q == 8) {
						omWE[q]         = omWE[7];
						omWE[q+N_OMINT] = omWE[7+N_OMINT];
						omSN[q]         = omSN[7];
						omSN[q+N_OMINT] = omSN[7+N_OMINT];
						omBT[q]         = omBT[7];
						omBT[q+N_OMINT] = omBT[7+N_OMINT];
					} else {
#endif

						// left-handed, central and right-handed
						// gradients for all directions
						REAL omxp =  gdata.om[q](i+1,j  ,k  ) - gdata.om[q](i  ,j  ,k  );
						REAL omx0 = (gdata.om[q](i+1,j  ,k  ) -
						             gdata.om[q](i-1,j  ,k  ))/2.;
						REAL omxm =  gdata.om[q](i  ,j  ,k  ) - gdata.om[q](i-1,j  ,k  );
						REAL omyp =  gdata.om[q](i  ,j+1,k  ) - gdata.om[q](i  ,j  ,k  );
						REAL omy0 = (gdata.om[q](i  ,j+1,k  ) - 
						             gdata.om[q](i  ,j-1,k  ))/2.;
						REAL omym =  gdata.om[q](i  ,j  ,k  ) - gdata.om[q](i  ,j-1,k  );
						REAL omzp =  gdata.om[q](i  ,j  ,k+1) - gdata.om[q](i  ,j  ,k  );
						REAL omz0 = (gdata.om[q](i  ,j  ,k+1) -
						             gdata.om[q](i  ,j  ,k-1))/2.;
						REAL omzm =  gdata.om[q](i  ,j  ,k  ) - gdata.om[q](i  ,j  ,k-1);

	  
						A(0,0,0) =  gdata.om[q](i,j,k);
	    
						// gradient in x-direction
						//	    A(1,0,0) = Limiter.compute(omxp,omx0,omxm);
						A(1,0,0) = Limiter->compute(omxp,omx0,omxm);
						// gradient in y-direction
						A(0,1,0) = Limiter->compute(omyp,omy0,omym);
						// gradient in z-direction
						A(0,0,1) = Limiter->compute(omzp,omz0,omzm);

						// Here N_OMINT signifies the values at the upper boundary
						// point values on lower (higher) x-boundary
						omWE[q] = A(0,0,0) - (0.5+shiftWE)*A(1,0,0);

						omWE[q+N_OMINT] = omWE[q] + A(1,0,0);
						domdx[q] = A(1,0,0);
						// point values on lower (higher) y-boundary
						omSN[q] = A(0,0,0) - (0.5+shiftSN)*A(0,1,0);

						omSN[q+N_OMINT] = omSN[q] + A(0,1,0);
						domdy[q] = A(0,1,0);
						// point values on lower (higher) z-boundary
						omBT[q] = A(0,0,0) - (0.5+shiftBT)*A(0,0,1);
						omBT[q+N_OMINT] = omBT[q] + A(0,0,1);
						domdz[q] = A(0,0,1);
					}
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
				}
#endif

				if(i>=-1 && j>=-1 && k>=-1){
					// Loop over lower (pos = 0) and upper boundary (pos = 1)
					for (int pos = 0; pos <= 1; ++pos){

						int dir = N_OMINT*pos;

// ---------------------------------------------------------------	      
//      Flux-computations:
//----------------------------------------------------------------

// ---------------------------------------------------------------	      
//      x-direction:
//----------------------------------------------------------------

						REAL rhoinv = 1./omWE[0+dir];
	    
						// Computation of pressure
						REAL psq,Bsq; // Square of momentum (mag ind.)
						if(N_OMINT >= 8){
							psq = (sqr(omWE[1+dir]) +
							       sqr(omWE[2+dir]) +
							       sqr(omWE[3+dir]))*sqr(omWE[0+dir]);
							Bsq = (sqr(omWE[4+dir]) +
							       sqr(omWE[5+dir]) +
							       sqr(omWE[6+dir]));
							if(thermal) {
								pres = (yama-1)*omWE[7+dir];
								omWE[7+dir] = TransEth2E(rhoinv,psq,Bsq,omWE[7+dir]);
							} else {
								pres = omWE[0+dir]*omWE[7+dir];
								omWE[7+dir] = TransT2E(rhoinv,psq,Bsq,omWE[7+dir]);
							}
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
#if(AUX_ENERGY == ENTROPY)
							if(N_OMINT==9) {
								// Transform to entropy
								omWE[8+dir] = omWE[8+dir]*pow(omWE[0+dir],1-Problem.gamma);
							}
#endif
#endif
						} else {
							//	      pres = pressure(omWE[0+dir]);
							pres = eos->pressure(gdata, Problem, omWE[0+dir], i-0.5+pos, j, k);
						}
	    
						fluxWE[0+dir] = omWE[0+dir]*omWE[1+dir];

						fluxWE[1+dir] = omWE[0+dir]*sqr(omWE[1+dir])
							+0.5*(-sqr(omWE[4+dir])
							      +sqr(omWE[5+dir])
							      +sqr(omWE[6+dir]))+pres;

						fluxWE[2+dir] = +omWE[0+dir]*omWE[1+dir]*omWE[2+dir]
							-omWE[4+dir]*omWE[5+dir];

#ifdef ANGULAR_MOMENTUM
						fluxWE[2+dir] *= gdata.h1(i-0.5+pos, j,k);
#endif
	    
						fluxWE[3+dir] = +omWE[0+dir]*omWE[1+dir]*omWE[3+dir]
							-omWE[4+dir]*omWE[6+dir];
	    
						if(N_OMINT >= 8){
							fluxWE[7+dir] = ((omWE[7+dir] + pres + 0.5*Bsq)*omWE[1+dir] -
							                 (omWE[1+dir]*omWE[4+dir] +
							                  omWE[2+dir]*omWE[5+dir] +
							                  omWE[3+dir]*omWE[6+dir])*omWE[4+dir]);
#ifdef CONDUCTION
							TPlus  = (Problem.gamma-1.)*om[7](i+pos  ,j,k)/om[0](i+pos  ,j,k);
							TMinus = (Problem.gamma-1.)*om[7](i+pos-1,j,k)/om[0](i+pos-1,j,k);
							fluxWE[7+dir] -= Problem.kappa*(TPlus - TMinus)*2.*hx[0];
#endif
						}

						if(N_OMINT == 9) {
							fluxWE[8+dir] = (omWE[8+dir]*omWE[1+dir]);
#if(AUX_ENERGY != ENTROPY)
#ifdef CONDUCTION
							if(N_OMINT == 9) {
								fluxWE[8+dir] -= Problem.kappa*(TPlus - TMinus)*2.*hx[0];
							}
#endif
#endif
						}
	    

// ---------------------------------------------------------------	      
//      y-direction:
//----------------------------------------------------------------

						rhoinv = 1./omSN[0+dir];
	    
						if(N_OMINT >= 8){
							psq = (sqr(omSN[1+dir]) +
							       sqr(omSN[2+dir]) +
							       sqr(omSN[3+dir]))*sqr(omSN[0+dir]);
							Bsq = (sqr(omSN[4+dir]) +
							       sqr(omSN[5+dir]) +
							       sqr(omSN[6+dir]));
							if(thermal) {
								pres = (yama-1)*omSN[7+dir];
								omSN[7+dir] = TransEth2E(rhoinv,psq,Bsq,omSN[7+dir]);
							} else {
								pres = omSN[0+dir]*omSN[7+dir];
								omSN[7+dir] = TransT2E(rhoinv,psq,Bsq,omSN[7+dir]);
							}
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
#if(AUX_ENERGY == ENTROPY)
							if(N_OMINT == 9) {
								omSN[8+dir] = omSN[8+dir]*pow(omSN[0+dir],1-Problem.gamma);
							}
#endif
#endif
						} else {
							//	      pres = pressure(omSN[0+dir]);
							pres = eos->pressure(gdata, Problem, omSN[0+dir], i, j-0.5+pos, k);
						}
		
						fluxSN[0+dir] = omSN[0+dir]*omSN[2+dir];
	    
						fluxSN[1+dir] = +omSN[0+dir]*omSN[1+dir]*omSN[2+dir]
							-omSN[4+dir]*omSN[5+dir];

						fluxSN[2+dir] = omSN[0+dir]*sqr(omSN[2+dir])
							+0.5*(+sqr(omSN[4+dir])
							      -sqr(omSN[5+dir])
							      +sqr(omSN[6+dir]))+pres;

#ifdef ANGULAR_MOMENTUM
						fluxSN[2+dir] *= gdata.h1(i,j-0.5+pos,k);
#endif
	      
						fluxSN[3+dir] = +omSN[0+dir]*omSN[2+dir]*omSN[3+dir]
							-omSN[5+dir]*omSN[6+dir];
	      
						if(N_OMINT >= 8){
							fluxSN[7+dir] = ((omSN[7+dir] + pres + 0.5*(Bsq))*omSN[2+dir]
							                 -(omSN[1+dir]*omSN[4+dir] +
							                   omSN[2+dir]*omSN[5+dir] +
							                   omSN[3+dir]*omSN[6+dir])*omSN[5+dir]);
#ifdef CONDUCTION
							TPlus  = (yama-1.)*om[7](i,j+pos  ,k)/om[0](i,j+pos  ,k);
							TMinus = (yama-1.)*om[7](i,j+pos-1,k)/om[0](i,j+pos-1,k);
							fluxSN[7+dir] -= Problem.kappa*(TPlus - TMinus)*2.*hx[1];
#endif
						}

						if(N_OMINT == 9) {
							fluxSN[8+dir] = (omSN[8+dir]*omSN[2+dir]);
#ifdef CONDUCTION
#if(AUX_ENERGY != ENTROPY)
							fluxSN[8+dir] -= Problem.kappa*(TPlus - TMinus)*2.*hx[1];
#endif
#endif
						}

// ---------------------------------------------------------------	      
//      z-direction:
//----------------------------------------------------------------

						rhoinv = 1./omBT[0+dir];

						if(N_OMINT >= 8){
							psq = (sqr(omBT[1+dir]) +
							       sqr(omBT[2+dir]) +
							       sqr(omBT[3+dir]))*sqr(omBT[0+dir]);
							Bsq = (sqr(omBT[4+dir]) +
							       sqr(omBT[5+dir]) +
							       sqr(omBT[6+dir]));
							if(thermal) {
								pres = (yama-1)*omBT[7+dir];
								omBT[7+dir] = TransEth2E(rhoinv,psq,Bsq,omBT[7+dir]);
							} else {
								pres = omBT[0+dir]*omBT[7+dir];
								omBT[7+dir] = TransT2E(rhoinv,psq,Bsq,omBT[7+dir]);
							}
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
#if(AUX_ENERGY == ENTROPY)
							if(N_OMINT == 9) {
								omBT[8+dir] = omBT[8+dir]*pow(omBT[0+dir],1.-Problem.gamma);
							}
#endif
#endif
						} else {
							//	      pres = pressure(omBT[0+dir]);
							pres = eos->pressure(gdata, Problem, omBT[0+dir], i, j, k-0.5+pos);
						}
	      
						fluxBT[0+dir] = omBT[0+dir]*omBT[3+dir];

						fluxBT[1+dir] = +omBT[0+dir]*omBT[1+dir]*omBT[3+dir]
							-omBT[4+dir]*omBT[6+dir];
	    
						fluxBT[2+dir] = +omBT[0+dir]*omBT[2+dir]*omBT[3+dir]
							-omBT[5+dir]*omBT[6+dir];

#ifdef ANGULAR_MOMENTUM
						fluxBT[2+dir] *= gdata.h1(i,j,k-0.5+pos);
#endif
	    
						fluxBT[3+dir] = omBT[0+dir]*sqr(omBT[3+dir])
							+0.5*(+sqr(omBT[4+dir])
							      +sqr(omBT[5+dir])
							      -sqr(omBT[6+dir]))+pres;
		
						if(N_OMINT >= 8){
							fluxBT[7+dir] = ((omBT[7+dir] + pres + 0.5*(Bsq))*omBT[3+dir]
							                 -(omBT[1+dir]*omBT[4+dir] +
							                   omBT[2+dir]*omBT[5+dir] +
							                   omBT[3+dir]*omBT[6+dir])*omBT[6+dir]);
#ifdef CONDUCTION
							TPlus  = (yama-1.)*om[7](i,j,k+pos  )/om[0](i,j,k+pos  );
							TMinus = (yama-1.)*om[7](i,j,k+pos-1)/om[0](i,j,k+pos-1);
							fluxBT[7+dir] -= Problem.kappa*(TPlus - TMinus)*2.*hx[2];

#endif	    
						}
						if(N_OMINT == 9) {
							fluxBT[8+dir] = (omBT[8+dir]*omBT[3+dir]);
#ifdef CONDUCTION
#if(AUX_ENERGY != ENTROPY)
							fluxBT[8+dir] -= Problem.kappa*(TPlus - TMinus)*2.*hx[2];
#endif
#endif
						}

// ---------------------------------------------------------------	      
//      End of pos-loop
//----------------------------------------------------------------     
					}

//----------------------------------------------------------------     
//      Computation of local electric field
//----------------------------------------------------------------     

					if(Problem.mag) {

						ExSB = -(omSB[2]*omSB[6] - omSB[3]*omSB[5]);
						ExNB = -(omNB[2]*omNB[6] - omNB[3]*omNB[5]);
						ExST = -(omST[2]*omST[6] - omST[3]*omST[5]);
						ExNT = -(omNT[2]*omNT[6] - omNT[3]*omNT[5]);
	    
						EyWB = -(omWB[3]*omWB[4] - omWB[1]*omWB[6]);
						EyEB = -(omEB[3]*omEB[4] - omEB[1]*omEB[6]);
						EyWT = -(omWT[3]*omWT[4] - omWT[1]*omWT[6]);
						EyET = -(omET[3]*omET[4] - omET[1]*omET[6]);

						EzWS = -(omWS[1]*omWS[5] - omWS[2]*omWS[4]);
						EzES = -(omES[1]*omES[5] - omES[2]*omES[4]);
						EzWN = -(omWN[1]*omWN[5] - omWN[2]*omWN[4]);
						EzEN = -(omEN[1]*omEN[5] - omEN[2]*omEN[4]);

					}

#ifdef PHYSDISS

// ---------------------------------------------------------------
// 	Computation of dissipative fluxes:
// ----------------------------------------------------------------

					REAL tau_W[4] = {0.}, tau_E[4] = {0.};
					REAL tau_S[4] = {0.}, tau_N[4] = {0.};
					REAL tau_B[4] = {0.}, tau_T[4] = {0.};
					REAL Bx, By, Bz;
					REAL dissEW, dissEE, dissES, dissEN, dissEB, dissET;
	  
// ---------------------------------------------------------------	      
//      x-direction:
//----------------------------------------------------------------
	
#ifdef GEOM  
#if GEOM != CARTESIAN
					dxinv[0] = 1./(gdata.dx[0]*gdata.h0(i-0.5,j,k));
					dxinv[1] = 1./(gdata.dx[1]*gdata.h1(i-0.5,j,k));
					dxinv[2] = 1./(gdata.dx[2]*gdata.h2(i-0.5,j,k));
#endif
#endif

					// perpendicular components of div_v
					double divPerp = (domdy[2]*dxinv[1] +
					                  domdz[3]*dxinv[2]); // dvydy + dvzdz


					// tau_xx
					tau_W[1] = (2.*(gdata.om[1](i  ,j,k) - gdata.om[1](i-1,j,k))*dxinv[0]-
					            divPerp)*twothirds;

					// tau_xy
					tau_W[2] = ((gdata.om[2](i  ,j,k) - gdata.om[2](i-1,j,k))*dxinv[0] +
					            domdy[1]*dxinv[1]);

					// tau_xz
					tau_W[3] = ((gdata.om[3](i  ,j,k) - gdata.om[3](i-1,j,k))*dxinv[0] +
					            domdz[1]*dxinv[2]);

#ifdef GEOM
#if GEOM != CARTESIAN
					sources->mod_Geom_Visc_WE(gdata, Problem, omWE, tau_W, i-0.5, j, k, 0);
#endif
#endif


					if(N_OMINT >= 8) {

						// Cell centered magnetic field
						Bx = (gdata.om[4](i,j,k) + gdata.om[4](i-1,j,k))*0.5;
						By = (gdata.om[5](i,j,k) + gdata.om[5](i,j-1,k))*0.5;
						Bz = (gdata.om[6](i,j,k) + gdata.om[6](i,j,k-1))*0.5;

						// Current density @ i-1/2
						REAL Jym = (domdz[4]*dxinv[2] - 
						            (0.5*(gdata.om[6](i,j,k  ) - gdata.om[6](i-1,j,k  ) +
						                  gdata.om[6](i,j,k-1) - gdata.om[6](i-1,j,k-1)))*dxinv[0]);
						REAL Jzm = ((0.5*(gdata.om[5](i,j  ,k) - gdata.om[5](i-1,j  ,k) +
						                  gdata.om[5](i,j-1,k) - gdata.om[5](i-1,j-1,k)))*dxinv[0] -
						            domdy[4]*dxinv[1]);

						dissEW = (Jym*Bz - Jzm*By);

					}


#ifdef GEOM
#if GEOM != CARTESIAN
					dxinv[0] = 1./(gdata.dx[0]*gdata.h0(i+0.5,j,k));
					dxinv[1] = 1./(gdata.dx[1]*gdata.h1(i+0.5,j,k));
					dxinv[2] = 1./(gdata.dx[2]*gdata.h2(i+0.5,j,k));
#endif
#endif

					// tau_xx
					tau_E[1] = (2.*(gdata.om[1](i+1,j,k) - gdata.om[1](i  ,j,k))*dxinv[0]-
					            divPerp)*twothirds;

					// tau_xy
					tau_E[2] = ((gdata.om[2](i+1,j,k) - gdata.om[2](i  ,j,k))*dxinv[0] +
					            domdy[1]*dxinv[1]);

					// tau_xz
					tau_E[3] = ((gdata.om[3](i+1,j,k) - gdata.om[3](i  ,j,k))*dxinv[0] +
					            domdz[1]*dxinv[2]);
	  

#ifdef GEOM
#if GEOM != CARTESIAN
					sources->mod_Geom_Visc_WE(gdata, Problem, omWE, tau_E, i+0.5, j, k, N_OMINT);
#endif
#endif	  

					if(N_OMINT >= 8) {

						// Current density @ i+1/2

						REAL Jyp = (domdz[4+N_OMINT]*dxinv[2] - 
						            (0.5*(gdata.om[6](i+1,j,k  ) - gdata.om[6](i,j,k  ) +
						                  gdata.om[6](i+1,j,k-1) - gdata.om[6](i,j,k-1)))*dxinv[0]);
						REAL Jzp = ((0.5*(gdata.om[5](i+1,j  ,k) - gdata.om[5](i,j  ,k) +
						                  gdata.om[5](i+1,j-1,k) - gdata.om[5](i,j-1,k)))*dxinv[0] -
						            domdy[4+N_OMINT]*dxinv[1]);

						dissEE = (Jyp*Bz - Jzp*By);

					}

// ---------------------------------------------------------------	      
//      y-direction:
//----------------------------------------------------------------

#ifdef GEOM
#if GEOM != CARTESIAN
					dxinv[0] = 1./(gdata.dx[0]*gdata.h0(i,j-0.5,k));
					dxinv[1] = 1./(gdata.dx[1]*gdata.h1(i,j-0.5,k));
					dxinv[2] = 1./(gdata.dx[2]*gdata.h2(i,j-0.5,k));
#endif
#endif

					divPerp = (domdx[1]*dxinv[0] +
					           domdz[3]*dxinv[2]); // dvxdx + dvzdz

					// tau_xy
					tau_S[1] = ((gdata.om[1](i,j  ,k) - gdata.om[1](i,j-1,k))*dxinv[1] +
					            domdx[2]*dxinv[0]);
	  
					// tau_yy
					tau_S[2] = (2.*(gdata.om[2](i,j  ,k) - gdata.om[2](i,j-1,k))*dxinv[1]-
					            divPerp)*twothirds;

					// tau_yz
					tau_S[3] = ((gdata.om[3](i,j  ,k) - gdata.om[3](i,j-1,k))*dxinv[1] +
					            domdz[2]*dxinv[2]);

#ifdef GEOM
#if GEOM != CARTESIAN
					sources->mod_Geom_Visc_SN(gdata, Problem, omSN, tau_S, i, j-0.5, k, 0);
#endif 
#endif

					if(N_OMINT >= 8) {

						// Current density @ j-1/2

						REAL Jzm = (domdx[5]*dxinv[0] -
						            (0.5*(gdata.om[4](i  ,j,k) - gdata.om[4](i  ,j-1,k) +
						                  gdata.om[4](i-1,j,k) - gdata.om[4](i-1,j-1,k)))*dxinv[1]);
						REAL Jxm = ((0.5*(gdata.om[6](i,j,k  ) - gdata.om[6](i,j-1,k  ) +
						                  gdata.om[6](i,j,k-1) - gdata.om[6](i,j-1,k-1)))*dxinv[1] -
						            domdz[5]*dxinv[2]);

						dissES = (Jzm*Bx - Jxm*Bz);

					}


#ifdef GEOM
#if GEOM != CARTESIAN
					dxinv[0] = 1./(gdata.dx[0]*gdata.h0(i,j+0.5,k));
					dxinv[1] = 1./(gdata.dx[1]*gdata.h1(i,j+0.5,k));
					dxinv[2] = 1./(gdata.dx[2]*gdata.h2(i,j+0.5,k));
#endif
#endif
					// tau_xy
					tau_N[1] = ((gdata.om[1](i,j+1,k) - gdata.om[1](i,j  ,k))*dxinv[1] +
					            domdx[2]*dxinv[0]);
	  
					// tau_yy
					tau_N[2] = (2.*(gdata.om[2](i,j+1,k) - gdata.om[2](i,j  ,k))*dxinv[1]-
					            divPerp)*twothirds;

					// tau_yz
					tau_T[3] = ((gdata.om[3](i,j+1,k) - gdata.om[3](i,j  ,k))*dxinv[1] +
					            domdz[2]*dxinv[2]);

#ifdef GEOM
#if GEOM != CARTESIAN
					sources->mod_Geom_Visc_SN(gdata, Problem, omSN, tau_N, i, j+0.5, k, N_OMINT);
#endif
#endif 

					if(N_OMINT >= 8) {

						// Current density @ j+1/2

						REAL Jzp = (domdx[5+N_OMINT]*dxinv[0] -
						            (0.5*(gdata.om[4](i  ,j+1,k) - gdata.om[4](i  ,j,k) +
						                  gdata.om[4](i-1,j+1,k) - gdata.om[4](i-1,j,k)))*dxinv[1]);

						REAL Jxp = ((0.5*(gdata.om[6](i,j+1,k  ) - gdata.om[6](i,j,k  ) +
						                  gdata.om[6](i,j+1,k-1) - gdata.om[6](i,j,k-1)))*dxinv[1] -
						            domdz[5+N_OMINT]*dxinv[2]);

						dissEN = (Jzp*Bx - Jxp*Bz);

					}

// ---------------------------------------------------------------	      
//      z-direction:
//----------------------------------------------------------------

#ifdef GEOM
#if GEOM != CARTESIAN
					dxinv[0] = 1./(gdata.dx[0]*gdata.h0(i,j,k-0.5));
					dxinv[1] = 1./(gdata.dx[1]*gdata.h1(i,j,k-0.5));
					dxinv[2] = 1./(gdata.dx[2]*gdata.h2(i,j,k-0.5));
#endif
#endif

					// perpendicular components of div_v
					divPerp = (domdx[1]*dxinv[0] +
					           domdx[2]*dxinv[1]);

					// tau_xz
					tau_B[1] = ((gdata.om[1](i,j,k  ) - gdata.om[1](i,j,k-1))*dxinv[2] +
					            domdx[3]*dxinv[0]);

					// tau_yz
					tau_B[2] = ((gdata.om[2](i,j,k  ) - gdata.om[2](i,j,k-1))*dxinv[2] +
					            domdy[3]*dxinv[1]);

					// tau_zz
					tau_B[3] = (2.*(gdata.om[3](i,j,k  ) - gdata.om[3](i,j,k-1))*dxinv[2]-
					            divPerp)*twothirds;


#ifdef GEOM
#if GEOM != CARTESIAN
					sources->mod_Geom_Visc_BT(gdata, Problem, omBT, tau_B, i, j, k-0.5, 0);
#endif
#endif

					if(N_OMINT >= 8) {

						// Current density @ k-1/2

						REAL Jxm = (domdy[6]*dxinv[1] -
						            (0.5*(gdata.om[5](i,j  ,k) - gdata.om[5](i,j  ,k-1) +
						                  gdata.om[5](i,j-1,k) - gdata.om[5](i,j-1,k-1)))*dxinv[2]);
						REAL Jym = ((0.5*(gdata.om[4](i  ,j,k) - gdata.om[4](i  ,j,k-1) +
						                  gdata.om[4](i-1,j,k) - gdata.om[4](i-1,j,k-1)))*dxinv[2] -
						            domdx[6]*dxinv[0]);

						dissEB = (Jxm*By - Jym*Bx);

					}

#ifdef GEOM
#if GEOM != CARTESIAN
					dxinv[0] = 1./(gdata.dx[0]*gdata.h0(i,j,k+0.5));
					dxinv[1] = 1./(gdata.dx[1]*gdata.h1(i,j,k+0.5));
					dxinv[2] = 1./(gdata.dx[2]*gdata.h2(i,j,k+0.5));
#endif
#endif
					// tau_xz
					tau_T[1] = ((gdata.om[1](i,j,k+1) - gdata.om[1](i,j,k  ))*dxinv[2] +
					            domdx[3]*dxinv[0]);

					// tau_yz
					tau_T[2] = ((gdata.om[2](i,j,k+1) - gdata.om[2](i,j,k  ))*dxinv[2] +
					            domdy[3]*dxinv[1]);
	  
					// tau_zz
					tau_T[3] = (2.*(gdata.om[3](i,j,k+1) - gdata.om[3](i,j,k  ))*dxinv[2]-
					            divPerp)*twothirds;

#ifdef GEOM
#if GEOM != CARTESIAN
					sources->mod_Geom_Visc_BT(gdata, Problem, omBT, tau_T, i, j, k+0.5, N_OMINT);
#endif
#endif

					if(N_OMINT >= 8) {

						// Current density @ k+1/2

						REAL Jxp = (domdy[6+N_OMINT]*dxinv[1] -
						            (0.5*(gdata.om[5](i,j  ,k+1) - gdata.om[5](i,j  ,k) +
						                  gdata.om[5](i,j-1,k+1) - gdata.om[5](i,j-1,k)))*dxinv[2]);

						REAL Jyp = ((0.5*(gdata.om[4](i  ,j,k+1) - gdata.om[4](i  ,j,k) +
						                  gdata.om[4](i-1,j,k+1) - gdata.om[4](i-1,j,k)))*dxinv[2] -
						            domdx[6+N_OMINT]*dxinv[0]);

						dissET = (Jxp*By - Jyp*Bx);

					}

#endif
// ---------------------------------------------------------------	      
//      End of dissipation
//      Now: conversion velocity -> momentum
//----------------------------------------------------------------

					for (int pos = 0; pos <= 1; ++pos){
						int dir = N_OMINT*pos;

						omWE[1+dir] = omWE[0+dir]*omWE[1+dir];
						omWE[2+dir] = omWE[0+dir]*omWE[2+dir];
						omWE[3+dir] = omWE[0+dir]*omWE[3+dir];

						omSN[1+dir] = omSN[0+dir]*omSN[1+dir];
						omSN[2+dir] = omSN[0+dir]*omSN[2+dir];
						omSN[3+dir] = omSN[0+dir]*omSN[3+dir];

						omBT[1+dir] = omBT[0+dir]*omBT[1+dir];
						omBT[2+dir] = omBT[0+dir]*omBT[2+dir];
						omBT[3+dir] = omBT[0+dir]*omBT[3+dir];

					}

// ---------------------------------------------------------------
// 	Estimate of maximum velocities:
// ----------------------------------------------------------------

// ---------------------------------------------------------------            
//      x-direction:
//----------------------------------------------------------------

					REAL rhoinv_p = 1./omWE[0];
					REAL rhoinv_m = 1./ommx[0];

					// Flow velocity
					REAL u_p = (omWE[1]*rhoinv_p);
					REAL u_m = (ommx[1]*rhoinv_m);

					REAL bsqr_p = (sqr(omWE[4]) +
					               sqr(omWE[5]) +
					               sqr(omWE[6]));
					REAL bsqr_m = (sqr(ommx[4]) +
					               sqr(ommx[5]) +
					               sqr(ommx[6]));

					// Alfven velocity (global)
					REAL va2_p  = bsqr_p*rhoinv_p;
					REAL va2_m  = bsqr_m*rhoinv_m;
					// Alfven velocity (directional)
					REAL vax2_p = sqr(omWE[4])*rhoinv_p;
					REAL vax2_m = sqr(ommx[4])*rhoinv_m;

					REAL pres_p;
					REAL pres_m;

					if(N_OMINT >= 8){
						REAL psq_p = (sqr(omWE[1]) +
						              sqr(omWE[2]) +
						              sqr(omWE[3]));
						REAL psq_m = (sqr(ommx[1]) +
						              sqr(ommx[2]) +
						              sqr(ommx[3]));
						pres_p = eos->pressure(rhoinv_p,psq_p,bsqr_p,omWE[7]);
						pres_m = eos->pressure(rhoinv_m,psq_m,bsqr_m,ommx[7]);
					} else {
						// 	    pres_p = pressure(omWE[0]);
						// 	    pres_m = pressure(ommx[0]);
						pres_p = eos->pressure(gdata, Problem, omWE[0], i-0.5, j, k);
						pres_m = eos->pressure(gdata, Problem, ommx[0], i-0.5, j, k);
					}
					// Sound speed
					REAL c2_p   = yama*pres_p*rhoinv_p;
					REAL c2_m   = yama*pres_m*rhoinv_m;
					// Fast mode speed 
					REAL vf_p = norm*sqrt(va2_p+c2_p + 
					                      sqrt(sqr(va2_p+c2_p)-4.*c2_p*vax2_p));
					REAL vf_m = norm*sqrt(va2_m+c2_m + 
					                      sqrt(sqr(va2_m+c2_m)-4.*c2_m*vax2_m));

					v_max_p[0] = std::max(std::max(vf_p+u_p,vf_m+u_m),0.);
					v_max_m[0] = std::max(std::max(vf_p-u_p,vf_m-u_m),0.);

					REAL vmax = std::max(v_max_p[0],v_max_m[0]);

					// Local computation of cfl number
					REAL cfl_loc = vmax/gdata.dx[0];
#ifdef GEOM
#if GEOM != CARTESIAN
					cfl_loc /= gdata.h0(i-0.5,j,k);
#endif
#endif
					cfl_lin = std::max(cfl_lin, cfl_loc);

#ifdef PHYSDISS
					cfl_eta_loc = 4.*Problem.nu(gdata, i-0.5,j,k)/(sqr(gdata.dx[0]*gdata.h0(i-0.5,j,k)));
					cfl_eta = std::max(cfl_eta, cfl_eta_loc);
#endif

#ifdef CONDUCTION
					REAL Temp_p = pres_p*rhoinv_p;
					REAL Temp_m = pres_m*rhoinv_m;

					REAL fac_max = std::max(Temp_p/omWE[7], Temp_m/ommx[7]);

					cfl_eta_loc = 2.*kappa*sqr(gdata.hx[0])/sqr(gdata.h0(i-0.5,j,k))*fac_max;
					cfl_eta = std::max(cfl_eta, cfl_eta_loc);
#endif

// ---------------------------------------------------------------            
//      y-direction:
//----------------------------------------------------------------	

					rhoinv_p = 1./omSN[0];
					rhoinv_m = 1./ommy[0](i);

					u_p = (omSN[2]*rhoinv_p);
					u_m = (ommy[2](i)*rhoinv_m);

					bsqr_p = (sqr(omSN[4]) +
					          sqr(omSN[5]) +
					          sqr(omSN[6]));
					bsqr_m = (sqr(ommy[4](i)) +
					          sqr(ommy[5](i)) +
					          sqr(ommy[6](i)));

					va2_p  = bsqr_p*rhoinv_p;
					va2_m  = bsqr_m*rhoinv_m;

					REAL vay2_p = sqr(omSN[5])*rhoinv_p;
					REAL vay2_m = sqr(ommy[5](i))*rhoinv_m;

					if(N_OMINT >= 8){
						REAL psq_p = (sqr(omSN[1]) +
						              sqr(omSN[2]) +
						              sqr(omSN[3]));
						REAL psq_m = (sqr(ommy[1](i)) +
						              sqr(ommy[2](i)) +
						              sqr(ommy[3](i)));
						pres_p = eos->pressure(rhoinv_p,psq_p,bsqr_p,omSN[7]);
						pres_m = eos->pressure(rhoinv_m,psq_m,bsqr_m,ommy[7](i));
					} else {
						//             pres_p = pressure(omSN[0]);
						//             pres_m = pressure(ommy[0](i));
						pres_p = eos->pressure(gdata, Problem, omSN[0], i, j-0.5, k);
						pres_m = eos->pressure(gdata, Problem, ommy[0](i), i, j-0.5, k);
					}


					c2_p   = yama*pres_p*rhoinv_p;
					c2_m   = yama*pres_m*rhoinv_m;

					vf_p = norm*sqrt(va2_p+c2_p+sqrt(sqr(va2_p+c2_p)-4.*c2_p*vay2_p));
					vf_m = norm*sqrt(va2_m+c2_m+sqrt(sqr(va2_m+c2_m)-4.*c2_m*vay2_m));

					v_max_p[1] = std::max(std::max(vf_p+u_p,vf_m+u_m),0.);
					v_max_m[1] = std::max(std::max(vf_p-u_p,vf_m-u_m),0.);

					vmax = std::max(v_max_p[1],v_max_m[1]);
					cfl_loc = vmax/gdata.dx[1];
#ifdef GEOM
#if GEOM != CARTESIAN
					cfl_loc /= gdata.h1(i,j-0.5,k);
#endif
#endif
					cfl_lin = std::max(cfl_lin, cfl_loc);
          
#ifdef PHYSDISS
					cfl_eta_loc = 4.*Problem.nu(gdata, i,j-0.5,k)/(sqr(gdata.dx[1]*gdata.h1(i,j-0.5,k)));
					cfl_eta = std::max(cfl_eta, cfl_eta_loc);
#endif
          
#ifdef CONDUCTION
					Temp_p = pres_p*rhoinv_p;
					Temp_m = pres_m*rhoinv_m;

					fac_max = std::max(Temp_p/omSN[7], Temp_m/ommy[7](i));
          
					cfl_eta_loc = 2.*kappa*sqr(gdata.hx[1])/sqr(gdata.h1(i,j-0.5,k))*fac_max;
					cfl_eta = std::max(cfl_eta, cfl_eta_loc);
#endif

// ---------------------------------------------------------------            
//      z-direction:
//----------------------------------------------------------------	
	
					rhoinv_p = 1./omBT[0];
					rhoinv_m = 1./ommz[0](i,j);

					u_p = (omBT[3]*rhoinv_p);
					u_m = (ommz[3](i,j)*rhoinv_m);

					bsqr_p = (sqr(omBT[4]) +
					          sqr(omBT[5]) +
					          sqr(omBT[6]));
					bsqr_m = (sqr(ommz[4](i,j)) +
					          sqr(ommz[5](i,j)) +
					          sqr(ommz[6](i,j)));

					va2_p  = bsqr_p*rhoinv_p;
					va2_m  = bsqr_m*rhoinv_m;

					REAL vaz2_p = sqr(omBT[6])*rhoinv_p;
					REAL vaz2_m = sqr(ommz[6](i,j))*rhoinv_m;
	  
					if(N_OMINT >= 8){
						REAL psq_p = (sqr(omBT[1]) +
						              sqr(omBT[2]) +
						              sqr(omBT[3]));
						REAL psq_m = (sqr(ommz[1](i,j)) +
						              sqr(ommz[2](i,j)) +
						              sqr(ommz[3](i,j)));
						pres_p = eos->pressure(rhoinv_p,psq_p,bsqr_p,omBT[7]);
						pres_m = eos->pressure(rhoinv_m,psq_m,bsqr_m,ommz[7](i,j));
					} else {
						//             pres_p = pressure(omBT[0]);
						//             pres_m = pressure(ommz[0](i,j));
						pres_p = eos->pressure(gdata, Problem, omBT[0], i, j, k-0.5);
						pres_m = eos->pressure(gdata, Problem, ommz[0](i,j), i, j, k-0.5);
					}


					c2_p   = yama*pres_p*rhoinv_p;
					c2_m   = yama*pres_m*rhoinv_m;

					vf_p = norm*sqrt(va2_p+c2_p+sqrt(sqr(va2_p+c2_p)-4.*c2_p*vaz2_p));
					vf_m = norm*sqrt(va2_m+c2_m+sqrt(sqr(va2_m+c2_m)-4.*c2_m*vaz2_m));

					v_max_p[2] = std::max(std::max(vf_p+u_p,vf_m+u_m),0.);
					v_max_m[2] = std::max(std::max(vf_p-u_p,vf_m-u_m),0.);

					vmax = std::max(v_max_p[2],v_max_m[2]);
	  
					cfl_loc  = vmax/gdata.dx[2];
#ifdef GEOM
#if GEOM != CARTESIAN
					cfl_loc /= gdata.h2(i,j,k-0.5);
#endif
#endif
					cfl_lin  = std::max(cfl_lin, cfl_loc);

#ifdef PHYSDISS
					cfl_eta_loc = 4.*Problem.nu(gdata, i,j,k-0.5)/(sqr(gdata.dx[2]*gdata.h2(i,j,k-0.5)));
					cfl_eta = std::max(cfl_eta, cfl_eta_loc);
#endif

#ifdef CONDUCTION
					Temp_p = pres_p*rhoinv_p;
					Temp_m = pres_m*rhoinv_m;

					fac_max = std::max(Temp_p/omBT[7], Temp_m/ommz[7](i,j));
	  
					cfl_eta_loc = 2.*kappa*sqr(gdata.hx[2])/sqr(gdata.h2(i,j,k-0.5))*fac_max;
					cfl_eta = std::max(cfl_eta, cfl_eta_loc);
#endif
  
// ---------------------------------------------------------------	      
//      End of dissipation
//      Now: conversion momentum -> angular momentum
//----------------------------------------------------------------

#ifdef ANGULAR_MOMENTUM
					for (int pos = 0; pos <= 1; ++pos){
						int dir = N_OMINT*pos;

						omWE[2+dir] *= gdata.h1(i-0.5+pos,j,k);
						omSN[2+dir] *= gdata.h1(i,j-0.5+pos,k);
						omBT[2+dir] *= gdata.h1(i,j,k-0.5+pos);
					}
#endif

// ---------------------------------------------------------------
//   Flux computation for all quantities
// ----------------------------------------------------------------

					REAL fac;
					REAL veps = 1.e-120;


					for (q = 0; q < N_OMINT; ++q){ 

						if(q < 4 || q > 6) {

							// 	    x-direction
	    
							fac = 2./(v_max_p[0]+v_max_m[0]+veps);

							HXpmWE[q](i,j,k) = (v_max_m[0]*fluxWE[q] + v_max_p[0]*flmx[q] -
							                    v_max_m[0]*v_max_p[0]*(omWE[q] - 
							                                           ommx[q]))*fac;

							// 	    y-direction
	    
							fac = 2./(v_max_p[1]+v_max_m[1]+veps);
	      
							HXpmSN[q](i,j,k) = (v_max_m[1]*fluxSN[q] + v_max_p[1]*flmy[q](i) -
							                    v_max_m[1]*v_max_p[1]*(omSN[q] -
							                                           ommy[q](i)))*fac;

	      
							// 	    z-direction
	    
							fac = 2./(v_max_p[2]+v_max_m[2]+veps);
	      
							HXpmBT[q](i,j,k) = (v_max_m[2]*fluxBT[q] + v_max_p[2]*flmz[q](i,j)-
							                    v_max_m[2]*v_max_p[2]*(omBT[q] - 
							                                           ommz[q](i,j)))*fac;
						}
					}

	  
// ---------------------------------------------------------------
//   Compute numerical emf directly from local emfs
// ----------------------------------------------------------------

					if(Problem.mag) {

						REAL BzBE_old(ommx[N_VAR+6]), BzBN_old(ommy[N_VAR+6](i));
						REAL BySE_old(ommx[N_VAR+5]), ByST_old(ommz[N_VAR+5](i,j)); 
						REAL BxWN_old(ommy[N_VAR+4](i)), BxWT_old(ommz[N_VAR+4](i,j));

						REAL ExST_old(flmz[4](i,j)), ExNB_old(flmy[4](i));
						REAL ExNT_old(flmz[N_OMINT+2](i,j-1));
						REAL EyWT_old(flmz[5](i,j)), EyEB_old(flmx[5]);
						REAL EyET_old(flmz[N_OMINT+1](i-1,j));
						REAL EzWN_old(flmy[6](i)), EzES_old(flmx[6]);
						REAL EzEN_old(flmy[N_OMINT](i-1));

						REAL vy_max_m_km(vym_old(i,j)), vy_max_p_km(vyp_old(i,j));
						REAL vz_max_m_jm(vzm_old(i)), vz_max_p_jm(vzp_old(i));
						REAL vx_max_m_km(vxm_old(i,j)), vx_max_p_km(vxp_old(i,j));
						REAL vz_max_m_im(vzm_old(i-1)), vz_max_p_im(vzp_old(i-1));
						REAL vx_max_m_jm(vxm_old(i,j-1)), vx_max_p_jm(vxp_old(i,j-1));
						REAL vy_max_m_im(vym_old(i-1,j)), vy_max_p_im(vyp_old(i-1,j));


						if(i >= 0 && j >= 0 && k >= 0) {

							v_cor_m[1] = std::max(v_max_m[1], vy_max_m_km);
							v_cor_p[1] = std::max(v_max_p[1], vy_max_p_km);
							v_cor_m[2] = std::max(v_max_m[2], vz_max_m_jm);
							v_cor_p[2] = std::max(v_max_p[2], vz_max_p_jm);
	      
							REAL facy = 1./(v_cor_p[1]+v_cor_m[1]+veps);
							REAL facz = 1./(v_cor_p[2]+v_cor_m[2]+veps);
							fac = facy*facz;
	      
							Ex(i,j,k) = ((v_cor_m[1]*v_cor_m[2]*ExSB +
							              v_cor_p[1]*v_cor_m[2]*ExNB_old +
							              v_cor_m[1]*v_cor_p[2]*ExST_old +
							              v_cor_p[1]*v_cor_p[2]*ExNT_old)*fac +
							             v_cor_m[1]*v_cor_p[1]*(omSB[6] - BzBN_old)*facy -
							             v_cor_m[2]*v_cor_p[2]*(omSB[5] - ByST_old)*facz);


							v_cor_m[0] = std::max(v_max_m[0], vx_max_m_km);
							v_cor_p[0] = std::max(v_max_p[0], vx_max_p_km);
							v_cor_m[2] = std::max(v_max_m[2], vz_max_m_im);
							v_cor_p[2] = std::max(v_max_p[2], vz_max_p_im);

							REAL facx = 1./(v_cor_p[0]+v_cor_m[0]+veps);
							facz = 1./(v_cor_p[2]+v_cor_m[2]+veps);
							fac = facx*facz;
	      
							Ey(i,j,k) = ((v_cor_m[0]*v_cor_m[2]*EyWB +
							              v_cor_p[0]*v_cor_m[2]*EyEB_old +
							              v_cor_m[0]*v_cor_p[2]*EyWT_old +
							              v_cor_p[0]*v_cor_p[2]*EyET_old)*fac +
							             v_cor_m[2]*v_cor_p[2]*(omWB[4] - BxWT_old)*facz -
							             v_cor_m[0]*v_cor_p[0]*(omWB[6] - BzBE_old)*facx);



							v_cor_m[0] = std::max(v_max_m[0], vx_max_m_jm);
							v_cor_p[0] = std::max(v_max_p[0], vx_max_p_jm);
							v_cor_m[1] = std::max(v_max_m[1], vy_max_m_im);
							v_cor_p[1] = std::max(v_max_p[1], vy_max_p_im);
	      
							facx = 1./(v_cor_p[0]+v_cor_m[0]+veps);
							facy = 1./(v_cor_p[1]+v_cor_m[1]+veps);
							fac = facx*facy;
	      
							Ez(i,j,k) = ((v_cor_m[0]*v_cor_m[1]*EzWS +
							              v_cor_p[0]*v_cor_m[1]*EzES_old +
							              v_cor_m[0]*v_cor_p[1]*EzWN_old +
							              v_cor_p[0]*v_cor_p[1]*EzEN_old)*fac +
							             v_cor_m[0]*v_cor_p[0]*(omWS[5] - BySE_old)*facx -
							             v_cor_m[1]*v_cor_p[1]*(omWS[4] - BxWN_old)*facy);


						}
					}

// ---------------------------------------------------------------
// 	  Compute dissipative fluxes and add them to numerical fluxes
// 	  Beware: 0.5 is contained in factor hx[dim]
// ---------------------------------------------------------------
	  
#ifdef PHYSDISS

					REAL mu_W = Problem.nu(gdata, i-0.5,j,k)*gdata.om[0](i,j,k);
					REAL mu_S = Problem.nu(gdata, i,j-0.5,k)*gdata.om[0](i,j,k);
					REAL mu_B = Problem.nu(gdata, i,j,k-0.5)*gdata.om[0](i,j,k);

					REAL mu_E = Problem.nu(gdata, i+0.5,j,k)*gdata.om[0](i,j,k);
					REAL mu_N = Problem.nu(gdata, i,j+0.5,k)*gdata.om[0](i,j,k);
					REAL mu_T = Problem.nu(gdata, i,j,k+0.5)*gdata.om[0](i,j,k);

					for (q = 1; q < 4; ++q) {

						HXpmWE[q](i  ,j,k) -= mu_W*tau_W[q];
						if(i < gdata.mx[0]+1){
							HXpmWE[q](i+1,j,k) -= mu_E*tau_E[q];
						}

						HXpmSN[q](i,j  ,k) -= mu_S*tau_S[q];
						if(j < gdata.mx[1]+1){
							HXpmSN[q](i,j+1,k) -= mu_N*tau_N[q];
						}

						HXpmBT[q](i,j,k  ) -= mu_B*tau_B[q];
						if(k < gdata.mx[2]+1){
							HXpmBT[q](i,j,k+1) -= mu_T*tau_T[q];
						}

					}

					if(N_OMINT >= 8) {

						REAL eta_W = Problem.eta(gdata, i-0.5,j,k);
						REAL eta_S = Problem.eta(gdata, i,j-0.5,k);
						REAL eta_B = Problem.eta(gdata, i,j,k-0.5);

						REAL eta_E = Problem.eta(gdata, i+0.5,j,k);
						REAL eta_N = Problem.eta(gdata, i,j+0.5,k);
						REAL eta_T = Problem.eta(gdata, i,j,k+0.5);

						// Influence of viscosity
						for(int q=1; q<4; ++q) {
							HXpmWE[7](i  ,j,k) += mu_W*gdata.om[1](i,j,k)*tau_W[q];
							if(i < gdata.mx[0]+1){
								HXpmWE[7](i+1,j,k) += mu_E*gdata.om[1](i,j,k)*tau_E[q];
							}
	      
							HXpmSN[7](i,j  ,k) += mu_S*gdata.om[2](i,j,k)*tau_S[q];
							if(j < gdata.mx[1]+1){
								HXpmSN[7](i,j+1,k) += mu_N*gdata.om[2](i,j,k)*tau_N[q];
							}

							HXpmBT[7](i,j,k  ) += mu_B*gdata.om[3](i,j,k)*tau_B[q];
							if(k < gdata.mx[2]+1){
								HXpmBT[7](i,j,k+1) += mu_T*gdata.om[3](i,j,k)*tau_T[q];
							}
						}

						// Influence of resistivity:
						HXpmWE[7](i  ,j,k) += eta_W*dissEW;
						if(i < gdata.mx[0]+1){
							HXpmWE[7](i+1,j,k) += eta_E*dissEE;
						}

						HXpmSN[7](i,j  ,k) += eta_S*dissES;
						if(j < gdata.mx[1]+1){
							HXpmSN[7](i,j+1,k) += eta_N*dissEN;
						}

						HXpmBT[7](i,j,k  ) += eta_B*dissEB;
						if(k < gdata.mx[2]+1){
							HXpmBT[7](i,j,k+1) += eta_T*dissET;
						}

					}
#endif
				} else {

					/*
					  Transformation for variables, which were outside the above loop:
					  velocity -> momentum
					  thermal energy -> overall energy
					  thermal energy -> entropy (if N_OMINT > 7 && ENTROPY)
					*/

					if(N_OMINT >= 8) {
						// x-direction
						REAL rhoinv = 1./omWE[0+N_OMINT];
						REAL psq = (sqr(omWE[1+N_OMINT]) +
						            sqr(omWE[2+N_OMINT]) +
						            sqr(omWE[3+N_OMINT]))*sqr(omWE[0+N_OMINT]);
						REAL Bsq = (sqr(omWE[4+N_OMINT]) +
						            sqr(omWE[5+N_OMINT]) +
						            sqr(omWE[6+N_OMINT]));
						if(thermal) {
							omWE[7+N_OMINT] = TransEth2E(rhoinv,psq,Bsq,omWE[7+N_OMINT]);
						} else {
							omWE[7+N_OMINT] = TransT2E(rhoinv,psq,Bsq,omWE[7+N_OMINT]);
						}

						// y-direction
						rhoinv = 1./omSN[0+N_OMINT];
						psq = (sqr(omSN[1+N_OMINT]) +
						       sqr(omSN[2+N_OMINT]) +
						       sqr(omSN[3+N_OMINT]))*sqr(omSN[0+N_OMINT]);
						Bsq = (sqr(omSN[4+N_OMINT]) +
						       sqr(omSN[5+N_OMINT]) +
						       sqr(omSN[6+N_OMINT]));
						if(thermal) {
							omSN[7+N_OMINT] = TransEth2E(rhoinv,psq,Bsq,omSN[7+N_OMINT]);
						} else {
							omSN[7+N_OMINT] = TransT2E(rhoinv,psq,Bsq,omSN[7+N_OMINT]);
						}

						// z-direction
						rhoinv = 1./omBT[0+N_OMINT];
						psq = (sqr(omBT[1+N_OMINT]) +
						       sqr(omBT[2+N_OMINT]) +
						       sqr(omBT[3+N_OMINT]))*sqr(omBT[0+N_OMINT]);
						Bsq = (sqr(omBT[4+N_OMINT]) +
						       sqr(omBT[5+N_OMINT]) +
						       sqr(omBT[6+N_OMINT]));
						if(thermal) {
							omBT[7+N_OMINT] = TransEth2E(rhoinv,psq,Bsq,omBT[7+N_OMINT]);
						} else {
							omBT[7+N_OMINT] = TransT2E(rhoinv,psq,Bsq,omBT[7+N_OMINT]);
						}

#if(AUX_ENERGY == ENTROPY)
						// Thermal energy -> Entropy
						if(N_OMINT >= 9) {
							omWE[8+N_OMINT] = omWE[8+N_OMINT]*pow(omWE[0+N_OMINT],
							                                      1-Problem.gamma);
							omSN[8+N_OMINT] = omSN[8+N_OMINT]*pow(omSN[0+N_OMINT],
							                                      1-Problem.gamma);
							omBT[8+N_OMINT] = omBT[8+N_OMINT]*pow(omBT[0+N_OMINT],
							                                      1-Problem.gamma);
						}
#endif
					}

					// Velocities -> momentum:
					// x-direction
					omWE[1+N_OMINT] = omWE[0+N_OMINT]*omWE[1+N_OMINT];
					omWE[2+N_OMINT] = omWE[0+N_OMINT]*omWE[2+N_OMINT];
					omWE[3+N_OMINT] = omWE[0+N_OMINT]*omWE[3+N_OMINT];
	  
					// y-direction
					omSN[1+N_OMINT] = omSN[0+N_OMINT]*omSN[1+N_OMINT];
					omSN[2+N_OMINT] = omSN[0+N_OMINT]*omSN[2+N_OMINT];
					omSN[3+N_OMINT] = omSN[0+N_OMINT]*omSN[3+N_OMINT];
	  
					// z-direction
					omBT[1+N_OMINT] = omBT[0+N_OMINT]*omBT[1+N_OMINT];
					omBT[2+N_OMINT] = omBT[0+N_OMINT]*omBT[2+N_OMINT];
					omBT[3+N_OMINT] = omBT[0+N_OMINT]*omBT[3+N_OMINT];

					// Momentum -> Angular momentum
#ifdef ANGULAR_MOMENTUM
					omWE[2+N_OMINT] *= gdata.h1(i+0.5,j,k);
					omSN[2+N_OMINT] *= gdata.h1(i,j+0.5,k);
					omBT[2+N_OMINT] *= gdata.h1(i,j,k+0.5);
#endif

				}
				for (q = 0; q < N_OMINT; ++q){
					ommx[q]      = omWE[q+N_OMINT];
					flmx[q]      = fluxWE[q+N_OMINT];
					ommy[q](i)   = omSN[q+N_OMINT];
					flmy[q](i)   = fluxSN[q+N_OMINT];
					ommz[q](i,j) = omBT[q+N_OMINT];
					flmz[q](i,j) = fluxBT[q+N_OMINT];
				}

				// Saving corner variables:
				ommx[N_VAR+5] = omES[5];
				ommx[N_VAR+6] = omEB[6];
				ommy[N_VAR+4](i)   = omWN[4];
				ommy[N_VAR+6](i)   = omNB[6];
				ommz[N_VAR+4](i,j) = omWT[4];
				ommz[N_VAR+5](i,j) = omST[5];

				flmx[5]      = EyEB;
				flmx[6]      = EzES;
				flmy[4](i)   = ExNB;
				flmy[6](i)   = EzWN;
				flmz[4](i,j) = ExST;
				flmz[5](i,j) = EyWT;
				flmz[6](i,j) = ExNT;
				flmz[N_OMINT](i,j) = EyET;
				flmy[5](i)   = EzEN;

				vzm_old(i) = v_max_m[2];
				vzp_old(i) = v_max_p[2];
				vxm_old(i,j) = v_max_m[0];
				vxp_old(i,j) = v_max_p[0];
				vym_old(i,j) = v_max_m[1];
				vyp_old(i,j) = v_max_p[1];

			}
		}
	}

	for(q=0; q<N_OMINT; ++q) {
		HXpmWE[q](-1,-1,-1) = 0.;
		HXpmSN[q](-1,-1,-1) = 0.;
		HXpmBT[q](-1,-1,-1) = 0.;
	}


#ifdef GEOM
#if GEOM != CARTESIAN
	// Modification of hyperbolic fluxes due to geometrical terms

	for (k = -1; k <= gdata.mx[2]+1; ++k){
		for (j = -1; j <= gdata.mx[1]+1; ++j){
			for (i = -1; i <= gdata.mx[0]+1; ++i){
				for (q = 0; q < N_OMINT; ++q){
					if(!(q >= 4 && q <=6)){
						REAL f_geom_x = gdata.h1(i-0.5,j,k)*gdata.h2(i-0.5,j,k);
						HXpmWE[q](i,j,k) *= f_geom_x;
					}
				}
			}
		}
	}


	for (k = -1; k <= gdata.mx[2]+1; ++k){
		for (j = -1; j <= gdata.mx[1]+1; ++j){
			for (i = -1; i <= gdata.mx[0]+1; ++i){
				for (q = 0; q < N_OMINT; ++q){
					if(!(q >= 4 && q <=6)){
						REAL f_geom_y = gdata.h0(i,j-0.5,k)*gdata.h2(i,j-0.5,k);
						HXpmSN[q](i,j,k) *= f_geom_y;
					}
				}
			}
		}
	}

	for (k = -1; k <= gdata.mx[2]+1; ++k){
		for (j = -1; j <= gdata.mx[1]+1; ++j){
			for (i = -1; i <= gdata.mx[0]+1; ++i){
				for (q = 0; q < N_OMINT; ++q){
					if(!(q >= 4 && q <=6)){
						REAL f_geom_z = gdata.h0(i,j,k-0.5)*gdata.h1(i,j,k-0.5);
						HXpmBT[q](i,j,k) *= f_geom_z;
					}
				}
			}
		}
	}
#endif
#endif	      
	  

	/* 
	   Compute boundary condition for fluxes if desired they have to be
	   specified in corresponding problem class:
	*/

	Problem.bc_Flux(gdata, HXpmWE, HXpmSN, HXpmBT);


// ---------------------------------------------------------------
//   Compute changes for hyperbolic variables
// ----------------------------------------------------------------

	for (k = 0; k <= gdata.mx[2]; ++k){
		for (j = 0; j <= gdata.mx[1]; ++j){
			for (i = 0; i <= gdata.mx[0]; ++i){
				for (q = 0; q < N_OMINT; ++q){

					REAL h_geom(1.);

#ifdef GEOM
#if GEOM != CARTESIAN
					h_geom /= (gdata.h0(i,j,k)*gdata.h1(i,j,k)*gdata.h2(i,j,k));
#endif 
#endif
					if(q < 4 || q > 6) {
						// 	    x-direction
						nom[q](i,j,k) += (HXpmWE[q](i+1,j,k) -
						                  HXpmWE[q](i  ,j,k))*gdata.hx[0]*h_geom;
						// 	    y-direction
						nom[q](i,j,k) += (HXpmSN[q](i,j+1,k) -
						                  HXpmSN[q](i,j  ,k))*gdata.hx[1]*h_geom;
						// 	    z-direction
						nom[q](i,j,k) += (HXpmBT[q](i,j,k+1) -
						                  HXpmBT[q](i,j,k  ))*gdata.hx[2]*h_geom;
					}
				}
			}
		}
	}


	// Freeing some memory
	for (q = 0; q < 4; ++q) {
		HXpmWE[q].resize(Index::set(0,0,0),Index::set(0,0,0));
		HXpmSN[q].resize(Index::set(0,0,0),Index::set(0,0,0));
		HXpmBT[q].resize(Index::set(0,0,0),Index::set(0,0,0));
	}

// ----------------------------------------------------------------
//   Flux computation for magnetic field
// ----------------------------------------------------------------

	if(Problem.mag) {
#ifdef PHYSDISS

		for (k = 0; k <= gdata.mx[2]+1; ++k){
			for (j = 0; j <= gdata.mx[1]+1; ++j){
				for (i = 0; i <= gdata.mx[0]+1; ++i){

					/*
					  Adding dissipation to the electric fields directly:
					  Ey(i,j,k) <-> Ey(i-1/2,j,k-1/2)
					*/
					REAL eta = Problem.eta(gdata, i,j-0.5,k-0.5);
#if GEOM > 1
					REAL f_geom = 1./(gdata.h1(i,j-0.5,k-0.5)*gdata.h2(i,j-0.5,k-0.5));

					Ex(i,j,k) += eta*((gdata.om[6](i,j  ,k-1)*gdata.h2(i,j  ,k-1.5) - 
					                   gdata.om[6](i,j-1,k-1)*
					                   gdata.h2(i,j-1,k-1.5))*2.*gdata.hx[1] -
					                  (gdata.om[5](i,j-1,k  )*gdata.h1(i,j-1.5,k  ) - 
					                   gdata.om[5](i,j-1,k-1)*
					                   gdata.h1(i,j-1.5,k-1))*2.*gdata.hx[2])*f_geom;
#else
					Ex(i,j,k) += eta*((gdata.om[6](i,j  ,k-1) - 
					                   gdata.om[6](i,j-1,k-1))*2.*gdata.hx[1] -
					                  (gdata.om[5](i,j-1,k  ) - 
					                   gdata.om[5](i,j-1,k-1))*2.*gdata.hx[2]);
#endif

					eta = Problem.eta(gdata, i-0.5,j,k-0.5);
#if GEOM > 1
					f_geom = 1./(gdata.h0(i-0.5,j,k-0.5)*gdata.h2(i-0.5,j,k-0.5));
	  
					Ey(i,j,k) += eta*((gdata.om[4](i-1,j,k  )*gdata.h0(i-1.5,j,k  ) - 
					                   gdata.om[4](i-1,j,k-1)*
					                   gdata.h0(i-1.5,j,k-1))*2.*gdata.hx[2] -
					                  (gdata.om[6](i  ,j,k-1)*gdata.h2(i  ,j,k-1.5) -
					                   gdata.om[6](i-1,j,k-1)*
					                   gdata.h2(i-1,j,k-1.5))*2.*gdata.hx[0])*f_geom;
#else
					Ey(i,j,k) += eta*((gdata.om[4](i-1,j,k  ) - 
					                   gdata.om[4](i-1,j,k-1))*2.*gdata.hx[2] -
					                  (gdata.om[6](i  ,j,k-1) -
					                   gdata.om[6](i-1,j,k-1))*2.*gdata.hx[0]);
#endif
	  
					eta = Problem.eta(gdata, i-0.5,j-0.5,k);
#if GEOM > 1
					f_geom = 1./(gdata.h0(i-0.5,j-0.5,k)*gdata.h1(i-0.5,j-0.5,k));

					Ez(i,j,k) += eta*((gdata.om[5](i  ,j-1,k)*gdata.h1(i  ,j-1.5,k) - 
					                   gdata.om[5](i-1,j-1,k)*
					                   gdata.h1(i-1,j-1.5,k))*2.*gdata.hx[0] -
					                  (gdata.om[4](i-1,j  ,k)*gdata.h0(i-1.5,j  ,k) -
					                   gdata.om[4](i-1,j-1,k)*
					                   gdata.h0(i-1.5,j-1,k))*2.*gdata.hx[1])*f_geom;
#else
					Ez(i,j,k) += eta*((gdata.om[5](i  ,j-1,k) - 
					                   gdata.om[5](i-1,j-1,k))*2.*gdata.hx[0] -
					                  (gdata.om[4](i-1,j  ,k) -
					                   gdata.om[4](i-1,j-1,k))*2.*gdata.hx[1]);
#endif

				}
			}
		}
#endif

    
		Problem.bc_emf(gdata, Ex, Ey, Ez);


		// If magnetic induction is treated directly -> numerical curl for noms

		if(IntegrateA) {
      
			for (k = 0; k <= gdata.mx[2]+1; ++k){
				for (j = 0; j <= gdata.mx[1]+1; ++j){
					for (i = 0; i <= gdata.mx[0]; ++i){
	    
						nom[4](i,j,k) = Ex(i,j,k);
	    
					}
				}
			}
      
			for (k = 0; k <= gdata.mx[2]+1; ++k){
				for (j = 0; j <= gdata.mx[1]; ++j){
					for (i = 0; i <= gdata.mx[0]+1; ++i){
	    
						nom[5](i,j,k) = Ey(i,j,k);
	    
					}
				}
			}
      
			for (k = 0; k <= gdata.mx[2]; ++k){
				for (j = 0; j <= gdata.mx[1]+1; ++j){
					for (i = 0; i <= gdata.mx[0]+1; ++i){

						nom[6](i,j,k) = Ez(i,j,k);
	    
					}
				}
			}
      
		} else {
      
      

#if GEOM > 1
			for (k = 0; k <= gdata.mx[2]; ++k){
				for (j = 0; j <= gdata.mx[1]; ++j){
					for (i = -1; i <= gdata.mx[0]; ++i){
	    
						REAL f_geom = 2./(gdata.h1(i+0.5,j,k)*gdata.h2(i+0.5,j,k));
						nom[4](i,j,k) = (+(gdata.h2(i+0.5,j+0.5,k)*Ez(i+1,j+1,k  ) - 
						                   gdata.h2(i+0.5,j-0.5,k)*Ez(i+1,j  ,k  ))*gdata.hx[1]
						                 -(gdata.h1(i+0.5,j,k+0.5)*Ey(i+1,j  ,k+1) -
						                   gdata.h1(i+0.5,j,k-0.5)*Ey(i+1,j  ,k  ))*gdata.hx[2])*f_geom;
	    
					}
				}
			}


			for (k = 0; k <= gdata.mx[2]; ++k){
				for (j = -1; j <= gdata.mx[1]; ++j){
					for (i = 0; i <= gdata.mx[0]; ++i){
	    
						REAL f_geom = 2./(gdata.h0(i,j+0.5,k)*gdata.h2(i,j+0.5,k));
						nom[5](i,j,k) = (+(gdata.h0(i,j+0.5,k+0.5)*Ex(i  ,j+1,k+1) - 
						                   gdata.h0(i,j+0.5,k-0.5)*Ex(i  ,j+1,k  ))*gdata.hx[2]
						                 -(gdata.h2(i+0.5,j+0.5,k)*Ez(i+1,j+1,k  ) - 
						                   gdata.h2(i-0.5,j+0.5,k)*Ez(i  ,j+1,k  ))*gdata.hx[0])*f_geom;
	    
	    
					}
				}
			}
			for (k = -1; k <= gdata.mx[2]; ++k){
				for (j = 0; j <= gdata.mx[1]; ++j){
					for (i = 0; i <= gdata.mx[0]; ++i){
	    
						REAL f_geom = 2./(gdata.h0(i,j,k+0.5)*gdata.h1(i,j,k+0.5));
						nom[6](i,j,k) = (+(gdata.h1(i+0.5,j,k+0.5)*Ey(i+1,j  ,k+1) - 
						                   gdata.h1(i-0.5,j,k+0.5)*Ey(i  ,j  ,k+1))*gdata.hx[0]
						                 -(gdata.h0(i,j+0.5,k+0.5)*Ex(i  ,j+1,k+1) - 
						                   gdata.h0(i,j-0.5,k+0.5)*Ex(i  ,j  ,k+1))*gdata.hx[1])*f_geom;
	    
					}
				}
			}
#else      
			for (k = 0; k <= gdata.mx[2]; ++k){
				for (j = 0; j <= gdata.mx[1]; ++j){
					for (i = -1; i <= gdata.mx[0]; ++i){
	    
						nom[4](i,j,k) = (+(Ez(i+1,j+1,k  ) - 
						                   Ez(i+1,j  ,k  ))*gdata.hx[1]
						                 -(Ey(i+1,j  ,k+1) -
						                   Ey(i+1,j  ,k  ))*gdata.hx[2])*2.;
	    
					}
				}
			}
    
			for (k = 0; k <= gdata.mx[2]; ++k){
				for (j = -1; j <= gdata.mx[1]; ++j){
					for (i = 0; i <= gdata.mx[0]; ++i){
	    
						nom[5](i,j,k) = (+(Ex(i  ,j+1,k+1) - 
						                   Ex(i  ,j+1,k  ))*gdata.hx[2]
						                 -(Ez(i+1,j+1,k  ) - 
						                   Ez(i  ,j+1,k  ))*gdata.hx[0])*2.;
	    
	    
					}
				}
			}
      
			for (k = -1; k <= gdata.mx[2]; ++k){
				for (j = 0; j <= gdata.mx[1]; ++j){
					for (i = 0; i <= gdata.mx[0]; ++i){
	    
						nom[6](i,j,k) = (+(Ey(i+1,j  ,k+1) - 
						                   Ey(i  ,j  ,k+1))*gdata.hx[0]
						                 -(Ex(i  ,j+1,k+1) - 
						                   Ex(i  ,j  ,k+1))*gdata.hx[1])*2.;
	    
					}
				}
			}
#endif
		}
	}

	for(q = 0; q<N_OMINT; ++q) {
		CheckNan(nom[q],q, 0, 1,"nom");
	}
  
// ----------------------------------------------------------------
//   Compute Courant number
// ----------------------------------------------------------------

#ifdef PHYSDISS
	if(cfl_eta > cfl_lin && n == RK_STEPS-1 && Problem.get_Info() && gdata.rank == 0) {
		cout << " Timestep given by dissipation " << endl;
	}
#else
	REAL cfl_eta(0.);
#endif
	REAL cfl = std::max(cfl_eta,cfl_lin);

#ifdef parallel
	double globalCFL = 0.;
	// Find maximum for CFL
	MPI_Barrier(gdata.comm3d);
	MPI_Reduce(&cfl, &globalCFL, 1, MPI_DOUBLE, MPI_MAX, 0, gdata.comm3d);
	MPI_Bcast(&globalCFL, 1, MPI_DOUBLE, 0, gdata.comm3d);
	cfl = globalCFL;
#endif
#ifdef DTCOMP_OLD
	cfl *= gdata.dt;
	if(gdata.rank==0) {
		if(n == RK_STEPS-1 && Problem.checkout(1)) {
			cout << "  cfl = " << cfl;
#ifdef PHYSDISS
			cout << " " << cfl_eta;
			cout << " " << comp_max_glob;
#endif
			cout << endl;
		}
	}
	if (cfl > 0.8) {
		cout << " *** cfl > cfl_max (=0.8) *** " << cfl << endl;
	}
#endif
#ifdef parallel
	MPI_Barrier(gdata.comm3d);
#endif

// ----------------------------------------------------------------
//   User defined source terms:
// ----------------------------------------------------------------

	Problem.src_User(gdata, nom, nom);

// ----------------------------------------------------------------
//   Geometrical source terms:
// ----------------------------------------------------------------

	sources->src_Geom(gdata, Problem, nom);

// ----------------------------------------------------------------
//   Other source terms (as, e.g., sources for thermal energy):
// ----------------------------------------------------------------

	sources->src_rhs(gdata, Problem, nom);

// ----------------------------------------------------------------
//   Check changes before applying them:
// ----------------------------------------------------------------

	for(q = 0; q<N_OMINT; ++q) {
		CheckNan(nom[q],q, 0, 2,"nom");
	}

// ----------------------------------------------------------------
//   Transformation to characteristic variables
// ----------------------------------------------------------------
 
	if(N_OMINT >= 8) {
		if(thermal) {
			Trafo->TransEth2E(gdata, gfunc, Problem);
		} else {
			Trafo->TransT2E(gdata, gfunc, Problem);
		}
	}

#ifdef ANGULAR_MOMENTUM
	Trafo->TransVel2AngMom(gdata, gfunc, Problem);
#else
	Trafo->TransVel2Momen(gdata, gfunc, Problem);
#endif

  //  TransMag2Pot(gdata, gfunc, Problem);

// ----------------------------------------------------------------
//   Determine domain to be integrated and apply changes:
// ----------------------------------------------------------------

//   int iend[3] = {gdata.mx[0], gdata.mx[1], gdata.mx[2]};

//   if(gfunc.get_bc_Type(1) == 1) {
//     iend[0] = gdata.mx[0]-1; // Shearing type boundaries
//   }


  // Domain: Bx: -1..mx[0], 0..mx[1], 0..mx[2];
  //          By: 0..mx[0], -1..mx[1], 0..mx[2];
  //          Bz: 0..mx[0], 0..mx[1], -1..mx[2];
  // Computing:  nom[4] = -((Ay(i+1,j,k+1) - Ay(i+1,j,k))*hx[2] -
  //                        (Az(i+1,j+1,k) - Az(i+1,j,k))*hx[1])
  //             nom[5] = -((Az(i+1,j+1,k) - Az(i,j+1,k))*hx[0] -
  //                        (Ax(i,j+1,k+1) - Ax(i,j+1,k))*hx[2])
  //             nom[6] = -((Ax(i,j+1,k+1) - Ax(i,j,k+1))*hx[1] -
  //                        (Ay(i+1,j,k+1) - Ay(i,j,k+1))*hx[0])
  // Therefore, we use the domain:
  //           Ax: 0..mx[0], 0..mx[1]+1, 0..mx[2]+1;
  //           Ay: 0..mx[0]+1, 0..mx[1], 0..mx[2]+1;
  //           Az: 0..mx[0]+1, 0..mx[1]+1, 0..mx[2];


	for (q = 0; q < N_OMINT; ++q){

		int qset(q);
		if(q >= 4 && q <= 6 && IntegrateA) {
			qset = N_OMINT+N_ADD+q-4;
		}

		int ibeg[3] = {0, 0, 0};
		int iend[3] = {gdata.mx[0], gdata.mx[1], gdata.mx[2]};

		if(q == 4) {
			if(IntegrateA) {
				iend[1] = gdata.mx[1]+1;
				iend[2] = gdata.mx[2]+1;
			} else {
				ibeg[0] = -1;
			}
		}
    
		if(q == 5) {
			if(IntegrateA) {
				iend[0] = gdata.mx[0]+1;
				iend[2] = gdata.mx[2]+1;
			} else {
				ibeg[1] = -1;
			}
		}

		if(q == 6) {
			if(IntegrateA) {
				iend[0] = gdata.mx[0]+1;
				iend[1] = gdata.mx[1]+1;
			} else {
				ibeg[2] = -1;
			}
		}

		if(gfunc.get_bc_Type(1) == 1) {
			iend[0] -= 1;
		}

#if (RK_STEPS == 3)
		if (n == 0) { // First Runge Kutta step
			for (k = ibeg[2]; k <= iend[2]; ++k){
				for (j = ibeg[1]; j <= iend[1]; ++j){
					for (i = ibeg[0]; i <= iend[0]; ++i){
						gdata.om[qset](i,j,k) -= gdata.dt*nom[q](i,j,k);
					}
				}
			}
		} else if (n == 1) { // Second Runge Kutta step
			for (k = ibeg[2]; k <= iend[2]; ++k) {
				for (j = ibeg[1]; j <= iend[1]; ++j) {
					for (i = ibeg[0]; i <= iend[0]; ++i) {
						gdata.om[qset](i,j,k) = 0.75*gdata.om[q+N_OMEGA](i,j,k)
							+0.25*gdata.om[qset](i,j,k)
							-0.25*gdata.dt*nom[q](i,j,k);
					}
				}
			}
		} else  if (n == 2) { // Third Runge Kutta step
			REAL twth   = 2./3.;
			REAL twthdt = twth*gdata.dt;

			for (k = ibeg[2]; k <= iend[2]; ++k){
				for (j = ibeg[1]; j <= iend[1]; ++j){
					for (i = ibeg[0]; i <= iend[0]; ++i){
						gdata.om[qset](i,j,k) = 1./3.*gdata.om[q+N_OMEGA](i,j,k)
							+twth*gdata.om[qset](i,j,k)
							-twthdt*nom[q](i,j,k);
					}
				}
			}
		}
#else
		if (n == 0) { // First Runge Kutta step
			for (k = ibeg[2]; k <= iend[2]; ++k){
				for (j = ibeg[1]; j <= iend[1]; ++j){
					for (i = ibeg[0]; i <= iend[0]; ++i){
						gdata.om[qset](i,j,k) -= gdata.dt*nom[q](i,j,k);
					}
				}
			}
		} else if (n == 1) { // Second Runge Kutta step
			for (k = ibeg[2]; k <= iend[2]; ++k) {
				for (j = ibeg[1]; j <= iend[1]; ++j) {
					for (i = ibeg[0]; i <= iend[0]; ++i) {
						gdata.om[qset](i,j,k) = 0.5*gdata.om[q+N_OMEGA](i,j,k)
							+0.5*gdata.om[qset](i,j,k)
							-0.5*gdata.dt*nom[q](i,j,k);
					}
				}
			}
		}
#endif
	}


	delete [] nom;

	gettimeofday(&tock, 0);
	cstep = clock();

// ----------------------------------------------------------------
//   Determine time at intermediate steps:
// ----------------------------------------------------------------

#if (RK_STEPS == 3) 
	// Time must come before boundary for time dependent boundaries
	if (n == 0) gdata.time += gdata.dt;      // First trial step with t* = t+dt
	if (n == 1) gdata.time -= 0.5*gdata.dt;  // Second trial step with t* = t+0.5dt
	if (n == 2) gdata.time += 0.5*gdata.dt;      // Last Runge-Kutta step - next time step
#else
	if (n == 0) gdata.time += gdata.dt;
#endif

// ----------------------------------------------------------------
//   Transformation to base variables and bcs
// ----------------------------------------------------------------

	gettimeofday(&tock2, 0);

#ifdef ANGULAR_MOMENTUM
	Trafo->TransAngMom2Vel(gdata, gfunc, Problem);
#else
	Trafo->TransMomen2Vel(gdata, gfunc, Problem);
#endif

	if(IntegrateA) {
		compute_B(gdata, gfunc, Problem);
	}

	for(int q=0; q<7; ++q) {
		gfunc.boundary(gdata, Problem, gdata.om[q],B,q);
	}


	if(N_OMINT >= 8) {

		if(thermal) {
			Trafo->TransE2Eth(gdata, gfunc, Problem, n, true);
		} else {
			Trafo->TransE2T(gdata, gfunc, Problem);
		}


		for(int q=7; q<N_OMINT; ++q) {
			gfunc.boundary(gdata, Problem, gdata.om[q],B,q);
		}
	}


// ----------------------------------------------------------------
//   Additional sources not included in nom.
// ----------------------------------------------------------------
  
	bool do_bcs[N_OMINT];
	for(int q=0; q<N_OMINT; ++q) {
		do_bcs[q] = false;
	}

	Problem.src_Add(gdata, do_bcs);

	for(int q=0; q<N_OMINT; ++q) {
		if(do_bcs[q]) {
			gfunc.boundary(gdata, Problem, gdata.om[q],B,q);
		}
	}

// ----------------------------------------------------------------
//   Performance estimate
// ----------------------------------------------------------------

	gettimeofday(&tock2, 0);

	timeval tstepOld = tstep;
	if(n == RK_STEPS-1) {
		gettimeofday(&tstep, 0);
	}
	cend = clock();

	REAL delt = ((tock.tv_sec + tock.tv_usec/1.e6) - 
	             (tick.tv_sec + tick.tv_usec/1.e6));

	REAL delt2 = ((tock2.tv_sec + tock2.tv_usec/1.e6) - 
	              (tock.tv_sec + tock.tv_usec/1.e6));

	REAL dtStep = ((tstep.tv_sec + tstep.tv_usec/1.e6) - 
	               (tstepOld.tv_sec + tstepOld.tv_usec/1.e6));

	if(gdata.rank == 0) {
		if(n == RK_STEPS-1 && Problem.get_Info() && Problem.checkout(5)) {
			cout << "------------------------------------------------------" << endl;
			cout << " Time needed for substeps: " << delt << " " << delt2 << endl;
			cout << " Full timestep:            " << dtStep << endl;
			//       cout << " CPU Cycle times: " << gdata.rank << " ";
			//       cout << (1.*(cstep-cstart))/CLOCKS_PER_SEC << " ";
			//       cout << (1.*(cend-cstep))/CLOCKS_PER_SEC << endl;
		}
	}

// ----------------------------------------------------------------
//   Test physical state
// ----------------------------------------------------------------

	phystest(gdata, gfunc, Problem, n);

	// Computing div B only at end of timestep
	if(Problem.mag && n == RK_STEPS-1){
		compute_divB(gdata, gfunc, Problem);
	}

	return cfl;
}

