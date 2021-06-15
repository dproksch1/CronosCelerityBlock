
#include "util.H"
#include<stdio.h>
#include <iomanip>

class semiLagra: public eulersolver {

public:
	void writesomething(Data &);
	semiLagra();
	semiLagra(Data &);
	REAL singlestep(Data &, gridFunc &, ProblemType &);
	void init(Data &, gridFunc &, ProblemType &);
	~semiLagra();

private:
	REAL norm(REAL);
	REAL vanLeer(REAL &,REAL &,REAL &);
	REAL DiLog (REAL);
	REAL findB(int, REAL);
	REAL Edot(int,REAL* ,REAL,int);
	REAL EccentricAnomaly(REAL,REAL);
	int getEnergy();
	//normalization
	REAL rho_norm, c_sound, time_norm;
	REAL gamma,Temp_0;
	//stars
	REAL M_Star1,R_Star1,T_Star1;
	REAL M_Star2,R_Star2,T_Star2;
	REAL Luminosity1, Luminosity2;
	REAL pos_comp1[DIM], pos_comp2[DIM];
	//wind
	REAL M_dot1, v_terminal1;
	REAL M_dot2, v_terminal2;
	//orbit
	REAL ROT_start;
	REAL semimajor, eccentricity, psi, apastron;
	REAL theta;
	//B-field
	REAL B_surface_Star1, rA_Star1, vrot_Star1, x_toro_Star1;
	REAL B_surface_Star2, rA_Star2, vrot_Star2, x_toro_Star2; 
	//fore-factors for energy densiy of radiation
	REAL fore_u_ph_1, fore_u_ph_2;
	// fore-factors for losses
	REAL fore_b_IC, fore_b_syn, fore_b_coul, fore_b_br, fore_b_coulP, fore_b_pp;
	REAL fore_acc;
	//energy bins: borders, center value, delta and correspond Klein-Nishina correction
	REAL Emax,Emin;
	REAL E[(int)(0.5*(N_OMINT_USER-6)+3)];
	REAL Ecent[(int)(0.5*(N_OMINT_USER-6)+3)];
	REAL deltaE[(int)(0.5*(N_OMINT_USER-6)+3)];
	REAL full_KN[(int)(0.5*(N_OMINT_USER-6)+3)];
	// global parameters
	REAL E_zero;
	REAL eps_T, T_elec,x_m,Zet;
	REAL diff_eff, inj_el, inj_pr;
	REAL M_p, M_e, M_pi;
	REAL AU_to_m;
// dummy output
	
			
};

semiLagra::~semiLagra(){
	
}

semiLagra::semiLagra(Data &gdata) : eulersolver(gdata) {

// NORMALIZATION
	// velocity normalization
	gamma = value((char*)"Adiabatic_exponent");
	Temp_0 = value((char*)"TempNorm");
	c_sound = sqrt(gamma*kB*Temp_0/m_H); // m s⁻¹
	rho_norm = 1e-8; // [kg m⁻³]
	time_norm = RSun/c_sound;
// STELLAR PARAMETERS	
	this->M_Star1 = value((char*)"Stellar_Mass1"); // [M_sol]
	this->R_Star1 = value((char*)"Stellar_Radius1"); // [RSun]
	this->T_Star1 = value((char*)"Teff1"); // [K]

	this->M_Star2 = value((char*)"Stellar_Mass2"); // [MSun]
	this->R_Star2 = value((char*)"Stellar_Radius2"); // [RSun]
	this->T_Star2 = value((char*)"Teff2"); // [K]

	// Some constants
	M_p = 938.272; // [MeV c²]  1,007276*AMU   proton mass in natural units
	M_pi = 139.570; // [MeV c²]  0.1499*AMU	pion mass in natural units
	M_e = m_e/echarge*sqr(c_light)/1e6;
	AU_to_m = AstronomicalUnit; 	

	Luminosity1 = 4*pi*stefan*sqr(R_Star1*RSun)*pow(T_Star1,4); // [J s⁻¹]
	Luminosity2 =4*pi*stefan*sqr(R_Star2*RSun)*pow(T_Star2,4); // [J s⁻¹]



// ORBITAL SETUP
	this->semimajor = value((char*)"semimajor"); // [RSun]
	this->eccentricity = value((char*)"eccentricity");
	ROT_start = value((char*)"Orbiting_after");
	apastron = value((char*)"apastron"); // 1 for ap- , 0 for periastron at t=0

	if (apastron ==1){psi = M_PI;}
	else if (apastron ==0){psi = 1e-5;}
		
	if (apastron ==1){
		pos_comp1[0] = -M_Star2/(M_Star1 + M_Star2)*semimajor*(1+eccentricity);
		pos_comp2[0] = +M_Star1/(M_Star1 + M_Star2)*semimajor*(1+eccentricity);
	} else if (apastron == 0){
		pos_comp1[0] = +M_Star2/(M_Star1 + M_Star2)*semimajor*(1-eccentricity);
		pos_comp2[0] = -M_Star1/(M_Star1 + M_Star2)*semimajor*(1-eccentricity);
	}
	pos_comp1[1] =  0;
	pos_comp1[2] =  0;
	pos_comp2[1] =  0;
	pos_comp2[2] =  0;	

// wind parameters
	this->M_dot1 = value((char*)"M_dot1"); // [MSun y⁻¹]
	this->v_terminal1 = value((char*)"v_terminal1"); //[m s⁻¹]
	
	this->M_dot2 = value((char*)"M_dot2"); //[MSun y⁻¹]	
	this->v_terminal2 = value((char*)"v_terminal2");// [m s⁻¹]		
	
// some assumptions on global parameters
	//some assumptions
	eps_T = 10; // [eV]
	B_surface_Star1 = 100; //G ... surface magentic field for star 1
	B_surface_Star2 = 100; //G ... surface magentic field for star 2
	diff_eff=1; // diffusion efficiency
	E_zero = 1; // [MeV]
	Zet = 1;
	inj_el=1e-5;
	inj_pr=1e-3;

		
// B-FIELD coefficients for findB(star,dist) following Usov and Melrose 1992
			
	// estimating Alfvén radius by setting ram pressure of stellar wind equal to magnetic stress component. 
	REAL xi_Star1= 0.1*sqr(B_surface_Star1)*sqr(R_Star1*RSun)
				/(M_dot1*MSun/year*v_terminal1);
	REAL xi_Star2= 0.1*sqr(B_surface_Star2)*sqr(R_Star2*RSun)
				/(M_dot2*MSun/year*v_terminal2);
					// Remark: 1 G² = 0.1 kg m⁻¹ s⁻²
	//Alfvén radius 
	if (xi_Star1>=1.){
		rA_Star1 = R_Star1*pow(xi_Star1,0.25); // [RSun]
	}else if (xi_Star1 <1.){
		rA_Star1 = R_Star1*(1.0+xi_Star1); // [RSun]
	}if (xi_Star2>=1.){
		rA_Star2 = R_Star2*pow(xi_Star2,0.25); // [RSun]
	}else if (xi_Star2 <1.){
		rA_Star2 = R_Star2*(1.0+xi_Star2); // [RSun]
	}

	// rotation velocity
	vrot_Star1 = 0.1*v_terminal1; // [m s⁻¹]
	vrot_Star2 = 0.1*v_terminal2; // [m s⁻¹]
	
	// condition for toroidal field
	x_toro_Star1 =R_Star1*(v_terminal1/vrot_Star1); // [RSun]
	x_toro_Star2 =R_Star2*(v_terminal2/vrot_Star2); // [RSun]
	
// LOSS-TERMS need the following factors
	// energy density of radiation	
	fore_u_ph_1 = (Luminosity1)/sqr(RSun)/(4*pi*c_light); 	// [J m⁻³]
	fore_u_ph_2 = (Luminosity2)/sqr(RSun)/(4*pi*c_light); 	// [J m⁻³]
	// IC cooling
	fore_b_IC = 4*sigma_Th/(3*sqr(m_e*c_light)*c_light)*MeV_to_Joule;
							// [MeV⁻¹ m³ s⁻¹ J⁻¹]
	// synchrotron emission
	fore_b_syn = 4*sigma_Th/(3*sqr(m_e*c_light)*c_light)*MeV_to_Joule;
							// [MeV⁻¹ m³ s⁻¹ J⁻¹]
	// bremstrahlung
	fore_b_br = 2/pi*fineStructure*sigma_Th*c_light;	// [m³ s⁻¹]
	// Coulomb losses electrons
	fore_b_coul = 55.725*c_light*sigma_Th*m_e*sqr(c_light)/MeV_to_Joule; // [m³ s⁻¹ MeV]
	// nucleon-nucleon-interaction
	fore_b_pp = 1.3*3*c_light*sigma_pp*M_pi/M_p; // [m³ s⁻¹]
	// Coulomb losses nucleons
	REAL CouLog= 20; // Coulomb logarithm
	T_elec = 1e8; //K
	x_m = 0.286*sqrt(T_elec/1e8);
	fore_b_coulP = 1.5*c_light*sigma_Th*m_e*sqr(c_light)/MeV_to_Joule*CouLog; // [m³s⁻¹ MeV]
		
// COMPUTE ENERGY BINS 	
	Emax = value((char*)"Emax");
	Emin = value((char*)"Emin");
	getEnergy();

	writesomething(gdata);

	

}


void semiLagra::init(Data &gdata, gridFunc & fu , ProblemType & fa ){
}



REAL semiLagra::singlestep(Data &gdata, gridFunc & gfunc , ProblemType & Problem){
	int bin;
#if (OMS_USER == TRUE)
// compute Kepler orbit
	if (gdata.time > ROT_start){
		
		psi = EccentricAnomaly(gdata.time-ROT_start,psi);
		theta = 2*atan(sqrt((1+eccentricity)/(1-eccentricity))*tan(psi/2));
		REAL distance = semimajor*(1-eccentricity*cos(psi));
// adjust stellar positions					
		pos_comp1[0] = M_Star2/(M_Star1 + M_Star2)*distance * cos(theta);
		pos_comp1[1] = M_Star2/(M_Star1 + M_Star2)*distance * sin(theta);
		pos_comp2[0] = -M_Star1/(M_Star1 + M_Star2)*distance * cos(theta);
		pos_comp2[1] = -M_Star1/(M_Star1 + M_Star2)*distance * sin(theta);
	}

//determine timestep for semiLagra
	REAL deltaT_global = gdata.dt*time_norm; //[s] // same as hydrocode
					//could do slope an several ...	

// LOOP OVER ALL CELLS		
	for(int k=0; k<=gdata.mx[2]; ++k) {
	for(int j=0; j<=gdata.mx[1]; ++j) {
	for(int i=0; i<=gdata.mx[0]; ++i) {

	//dummy output
		REAL out2=0;
		REAL out1=0;
		REAL out3=0;
			
		//int direc=0;		
		//int neigh[3];
		//neigh[0]=0;neigh[1]=0;neigh[2]=0;

	// get distance to stars
		REAL r_sph1 = sqrt(sqr(gdata.getCen_x(i) - pos_comp1[0]) +
			 	   sqr(gdata.getCen_y(j) - pos_comp1[1]) +
		                   sqr(gdata.getCen_z(k) - pos_comp1[2]));
		REAL r_sph2 = sqrt(sqr(gdata.getCen_x(i) - pos_comp2[0]) +
     	          	           sqr(gdata.getCen_y(j) - pos_comp2[1]) +
     		       	           sqr(gdata.getCen_z(k) - pos_comp2[2]));

	//NO ACC at edges:
		int disable = 0;
		// if (gdata.getCen_x(i)>value((char*)"xe")-2*gdata.dx[0] || gdata.getCen_x(i)<value((char*)"xb")+2*gdata.dx[0]){
		// 	disable = 1;
		// } else if (gdata.getCen_y(j)>value((char*)"ye")-2*gdata.dx[1] || gdata.getCen_y(j)<value((char*)"yb")+2*gdata.dx[1]){
		// 	disable = 1;
		// } else if (gdata.getCen_z(k)>value((char*)"ze")-2*gdata.dx[2]|| gdata.getCen_z(k)<value((char*)"zb")+2*gdata.dx[2]) {
		// 	disable = 1;
		// }

		if (gdata.getCen_x(i)>gdata.global_xe[0]-2*gdata.dx[0] ||
			gdata.getCen_x(i)<gdata.global_xb[0]+2*gdata.dx[0]){
			disable = 1;
		}
		if (gdata.getCen_y(j)>gdata.global_xe[1]-2*gdata.dx[1] ||
			gdata.getCen_y(j)<gdata.global_xb[1]+2*gdata.dx[1]){
			disable = 1;
		}
		if (gdata.getCen_z(k)>gdata.global_xe[2]-2*gdata.dx[2] ||
			gdata.getCen_z(k)<gdata.global_xb[2]+2*gdata.dx[2]) {
			disable = 1;
		}
		

	// Do acceleration loop only if we 1 cell above stellar surface
		if( r_sph1 >= R_Star1+4*gdata.dx[0] 
			&& r_sph2 >= R_Star2+4*gdata.dx[0]){
 
		// compute divergence of velocity field
			REAL dv_xdx=(gdata.om[q_sx](i+1,j,k)-gdata.om[q_sx](i-1,j,k))
					/(gdata.getCen_x(i+1)-gdata.getCen_x(i-1)); 
			REAL dv_ydy=(gdata.om[q_sy](i,j+1,k)-gdata.om[q_sy](i,j-1,k))
					/(gdata.getCen_y(j+1)-gdata.getCen_y(j-1));
			REAL dv_zdz=(gdata.om[q_sz](i,j,k+1)-gdata.om[q_sz](i,j,k-1))
					/(gdata.getCen_z(k+1)-gdata.getCen_z(k-1));
			REAL nablaV = dv_xdx + dv_ydy + dv_zdz; // [c_sound RSun⁻¹]

		
		// retrieve number density of particles in wind
			REAL N_H = gdata.om[q_rho](i,j,k)*rho_norm/m_H; //  [m⁻³]
		// get energy density of radiaton	
			REAL u_ph = fore_u_ph_1/sqr(r_sph1)+fore_u_ph_2/sqr(r_sph2); // [J m⁻³]
		// get IC-loss rate in Thomson-limit
			REAL b_IC = fore_b_IC*u_ph; 	// [MeV⁻¹ s⁻¹]
		// determine B-fiels	
			REAL B1 = findB(1,r_sph1); // [G] 
			REAL B2 = findB(2,r_sph2); // [G]
			REAL B_total = 1e-4*std::max(B1,B2); // [T]

			 // [T]
			//double E_Bohm_old = 16.0e15*3/1e6*1e-4*std::max(B1,B2);
			//double E_Bohm_new = 16.0e15*3/1e6*B_total;

	/*if ( abs(gdata.getCen_z(k))<0.1 && abs(gdata.getCen_y(j))<0.1){// && gdata.getCen_x(i)>-1000 && gdata.getCen_x(i)<-250 ){
			cout << "i: " << gdata.getCen_x(i) <<","<<gdata.getCen_y(j)<<","<<gdata.getCen_z(k) << "  rho: " << gdata.om[q_rho](i,j,k) << ", B: " << 1e-4*std::max(B1,B2) << " and new: " << B_total << endl;
cout << " E_Bohm: " << E_Bohm_old << " and new: " << E_Bohm_new << " MeV" << endl;
		}*/
				
				
		// determine energy density of magnetic fields
			REAL u_B = 0.1*sqr(B_total*1e4)/(8*pi); 
				//[G²]= [g cB_surface_Star1m⁻¹ s⁻²] = [erg cm⁻³] = 0.1 [J m⁻³]
		// synchrotron losses				
			REAL b_syn = fore_b_syn*u_B; // [MeV⁻¹ s⁻¹]
		// bremsstrahlung losses
			REAL b_br = fore_b_br*N_H; // [s⁻¹]
		// Coulomb losses electrons
			REAL b_coul = fore_b_coul*N_H; // [MeV s⁻¹]
		// Nucleon-Nucleon interaction losses
			REAL b_pp = fore_b_pp*N_H;  // [s⁻¹]
		// Coulomb losses protons
			REAL b_coulP = fore_b_coulP*sqr(Zet)*N_H; // [MeV s⁻¹]

		//ACCELERATION CELL
			REAL Edot_DSA = 0;  	// diffusive shock acceleraton
			//REAL fore_tau = 0; 
			REAL T_0 =0 ;		// diffusive escape timescale
			REAL E_Bohm = 0;	// Bohm_limit at energy
			REAL V_shock = 0;	// shock velocity
			REAL c_r=0;
		// identify acceleration cells by
			int acccell=0;
			if (nablaV<2 ){
				if ( gdata.om[q_Eges](i,j,k)<20){
					if (gdata.om[q_Eges](i+1,j,k)>20){acccell=1;}
					if (gdata.om[q_Eges](i-1,j,k)>20){acccell=1;}
					if (gdata.om[q_Eges](i,j+1,k)>20){acccell=1;}
					if (gdata.om[q_Eges](i,j-1,k)>20){acccell=1;}
					if (gdata.om[q_Eges](i,j,k+1)>20){acccell=1;}
					if (gdata.om[q_Eges](i,j,k-1)>20){acccell=1;}
				}
			}
			if (acccell==1 && disable ==0){
				out2 = 1;
				//direc=1;
			// compute gradient of temperature field
					REAL dTemp_dx=(gdata.om[q_Eges](i+3,j,k)-gdata.om[q_Eges](i-3,j,k))
							/(gdata.getCen_x(i+3)-gdata.getCen_x(i-3));
					REAL dTemp_dy= (gdata.om[q_Eges](i,j+3,k)-gdata.om[q_Eges](i,j-3,k))
							/(gdata.getCen_y(j+3)-gdata.getCen_y(j-3));
					REAL dTemp_dz= (gdata.om[q_Eges](i,j,k+3)-gdata.om[q_Eges](i,j,k-3))
							/(gdata.getCen_z(k+3)-gdata.getCen_z(k-3));
				//and its norm
					REAL dTempnorm = sqrt(sqr(dTemp_dx)+sqr(dTemp_dy)+sqr(dTemp_dz));
			//approximate compression ratio
				//identify neighbor cells along density gradient
				REAL Tempx = norm(dTemp_dx);
				REAL Tempy = norm(dTemp_dy);
				REAL Tempz = norm(dTemp_dz);
				/*if (Tempx>=Tempy && Tempx>=Tempz){ // neighbors along x
					neigh[0]=1; 
				}else if (Tempy>=Tempx && Tempy>=Tempz){
					neigh[1]=1; 
				}else if (Tempz>=Tempx && Tempz>=Tempy){
					neigh[2]=1; 
				}*/
				// get density of neighbor cells
				/*int thick=3;
				REAL rho_neigh1 = 
				gdata.om[q_rho](i+thick*neigh[0],j+thick*neigh[1],k+thick*neigh[2]);
				REAL rho_neigh2 = 
				gdata.om[q_rho](i-thick*neigh[0],j-thick*neigh[1],k-thick*neigh[2]);
				if (rho_neigh2 > rho_neigh1){direc=2;}
				REAL rho_max = std::max(rho_neigh1,rho_neigh2);
				REAL rho_min = std::min(rho_neigh1,rho_neigh2);
				if (gdata.om[q_rho](i,j,k)<rho_min){rho_min=gdata.om[q_rho](i,j,k);}*/

				REAL boint[3];
				boint[0]= gdata.getCen_x(i) + 0.5*gdata.dx[0];
				boint[1]= gdata.getCen_y(j) + 0.5*gdata.dx[1];
				boint[2]= gdata.getCen_z(k) + 0.5*gdata.dx[2];

				REAL stretch = 3*gdata.dx[0];
				REAL coint[3];
				coint[0]= boint[0] + stretch*dTemp_dx/dTempnorm;
				coint[1]= boint[1] + stretch*dTemp_dy/dTempnorm;
				coint[2]= boint[2] + stretch*dTemp_dz/dTempnorm;
			// corner points and weighting by distance
				REAL oned[3];
				oned[0] = (coint[0]-gdata.xb[0])/gdata.dx[0];
				oned[1] = (coint[1]-gdata.xb[1])/gdata.dx[1];
				oned[2] = (coint[2]-gdata.xb[2])/gdata.dx[2];
				int one[3];
				one[0] = (coint[0]-gdata.xb[0])/gdata.dx[0];
				one[1] = (coint[1]-gdata.xb[1])/gdata.dx[1];
				one[2] = (coint[2]-gdata.xb[2])/gdata.dx[2];
				if (oned[0]<=0){
					if (oned[0]>=-1){one[0]=-1;}
					else if (oned[0] >= -2){one[0]=-2;}
					else if (oned[0] >= -3){one[0]=-3;}
					else {cout << "ERRRRROR: " << oned[0] << " " << one[0] << endl;}
				}
				if (oned[1]<=0){
					if (oned[1]>=-1){one[1]=-1;}
					else if (oned[1] >= -2){one[1]=-2;}
					else if (oned[1] >= -3){one[1]=-3;}
					else {cout << "ERRRRROR: " << oned[1] << " " << one[1] << endl;}
				}
				if (oned[2]<=0){
					if (oned[2]>=-1){one[2]=-1;}
					else if (oned[2] >= -2){one[2]=-2;}
					else if (oned[2] >= -3){one[2]=-3;}
					else {cout << "ERRRRROR: " << oned[2] << " " << one[2] << endl;}
				}


				REAL rho_max = 0;
				REAL distsum = 0;
				REAL dist1 = gdata.dx[0]*sqrt(sqr(oned[0]-(one[0]+0.5))
						       +sqr(oned[1]-(one[1]+0.5))+
						        sqr(oned[2]-(one[2]+0.5)));
				if (dist1 < gdata.dx[0]){
					distsum += (1-dist1/gdata.dx[0]);//cout << (1-dist1/gdata.dx[0]) << " ";
					rho_max += (1-dist1/gdata.dx[0])*gdata.om[q_rho](one[0],one[1],one[2]);
				}
				REAL dist2 = gdata.dx[0]*sqrt(sqr(oned[0]-(one[0]+1.5))
						       +sqr(oned[1]-(one[1]+0.5))+
							sqr(oned[2]-(one[2]+0.5)));
				if (dist2 < gdata.dx[0]){
					distsum += (1-dist2/gdata.dx[0]);//cout << (1-dist2/gdata.dx[0]) << " ";
					rho_max += (1-dist2/gdata.dx[0])*gdata.om[q_rho](one[0]+1,one[1],one[2]);
				}
				REAL dist3 = gdata.dx[0]*sqrt(sqr(oned[0]-(one[0]-0.5))
						       +sqr(oned[1]-(one[1]+0.5))+
						        sqr(oned[2]-(one[2]+0.5)));
				if (dist3 < gdata.dx[0]){
					distsum += (1-dist3/gdata.dx[0]);//cout << (1-dist3/gdata.dx[0]) << " ";
					rho_max += (1-dist3/gdata.dx[0])*gdata.om[q_rho](one[0]-1,one[1],one[2]);
				}
				REAL dist4 = gdata.dx[0]*sqrt(sqr(oned[0]-(one[0]+0.5))
						       +sqr(oned[1]-(one[1]+1.5))+
							sqr(oned[2]-(one[2]+0.5)));
				if (dist4 < gdata.dx[0]){//cout << (1-dist4/gdata.dx[0]) << " ";
					distsum += (1-dist4/gdata.dx[0]);
					rho_max += (1-dist4/gdata.dx[0])*gdata.om[q_rho](one[0],one[1]+1,one[2]);
				}
				REAL dist5 = gdata.dx[0]*sqrt(sqr(oned[0]-(one[0]+0.5))
						       +sqr(oned[1]-(one[1]-0.5))+
							sqr(oned[2]-(one[2]+0.5)));
				if (dist5 < gdata.dx[0]){
					distsum += (1-dist5/gdata.dx[0]);//cout << (1-dist5/gdata.dx[0]) << " ";
					rho_max += (1-dist5/gdata.dx[0])*gdata.om[q_rho](one[0],one[1]-1,one[2]);
				}
				REAL dist6 = gdata.dx[0]*sqrt(sqr(oned[0]-(one[0]+0.5))
						       +sqr(oned[1]-(one[1]+0.5))+
							sqr(oned[2]-(one[2]+1.5)));
				if (dist6 < gdata.dx[0]){
					distsum += (1-dist6/gdata.dx[0]);//cout << (1-dist6/gdata.dx[0]) << " ";
					rho_max += (1-dist6/gdata.dx[0])*gdata.om[q_rho](one[0],one[1],one[2]+1);
				}
				REAL dist7 = gdata.dx[0]*sqrt(sqr(oned[0]-(one[0]+0.5))
						       +sqr(oned[1]-(one[1]+0.5))+
							sqr(oned[2]-(one[2]-0.5)));
				if (dist7 < gdata.dx[0]){
					distsum += (1-dist7/gdata.dx[0]);//cout << (1-dist7/gdata.dx[0]) << " ";
					rho_max += (1-dist7/gdata.dx[0])*gdata.om[q_rho](one[0],one[1],one[2]-1);
				}
				REAL rho_min = gdata.om[q_rho](i,j,k);
		
 		
				//cout << "direc: " << direc << endl;
				// determine compression
				c_r = rho_max/rho_min;
				//if (c_r>1.5){
				
				// approximate shock-velocity as v-component along the gradient of rho
					V_shock = (gdata.om[q_sx](i,j,k)*dTemp_dx 
					 	 	 + gdata.om[q_sy](i,j,k)*dTemp_dy 
					  	 	 + gdata.om[q_sz](i,j,k)*dTemp_dz)/dTempnorm;
					if(dTempnorm==0){V_shock=0;} //just in case
					V_shock *= c_sound; // [m s⁻¹]
				//}	
/*if (V_shock/1000 < 2232.52 && V_shock/1000 > 2232.50){
	cout << "INFO: " << endl;
	cout << i << " " << j << " " << k << endl;
	cout << "xb: " << gdata.xb[0] << " " << gdata.xb[1] << " " << gdata.xb[2] << endl;
	cout << "mb: " << gdata.mx[0] << " " << gdata.mx[1] << " " << gdata.mx[2] << endl;
	cout << "boint: " << boint[0] << " " << boint[1] << " " << boint[2] << endl;
	cout << "coint: " << coint[0] << " " << coint[1] << " " << coint[2] << endl;
	cout << "one: " << one[0] << " " << one[1] << " " << one[2] << endl;
	cout << "oned: " << oned[0] << " " << oned[1] << " " << oned[2] << endl;
	cout  << dist1 << " " ;
	cout  << dist2 << " " ;
	cout  << dist3 << " " ;
	cout  << dist4 << " " ;
	cout  << dist5 << " " ;
	cout  << dist6 << " " ;
	cout  << dist7 << " " ;
cout << endl ;
cout  << rho_max << endl
<< endl;


}*/
			// diffusion along shock
				REAL kappa_a = value((char*)"DiffusionConstant"); // [m² s⁻¹]
				// diffusion perpendicular to shock
				//REAL kappa_a = kappa_d*4;  // [m² s⁻¹]
				// convective flow
				//REAL pfr = 0.5;
				//REAL V_cone = pfr*V_shock; // [m s⁻¹]
				//diffusive timescale = convective timescale at
				//REAL r_0 = 4*kappa_d/V_cone; // [m]
				// diffusive escape timescale
				//T_0 = sqr(r_0)/4/kappa_d; // [s]
				
				// Bohm_limit at energy
				E_Bohm = kappa_a*B_total*3/1e6; // [MeV]

				// diffusive shock acceleraton
				if (V_shock > 0 ){
					T_0 = c_r*kappa_a/sqr(V_shock);
					Edot_DSA  = (c_r-1)/(3*c_r*kappa_a) * sqr(V_shock); // [s⁻¹]
				}
			// alternative
				// get energy-independant acceleration rate
				/*Edot_DSA = (c_r-1)/c_r/(c_r+1)*B_total/diff_eff
						* sqr(V_shock)/1e6; // [MeV s⁻¹]
				//get escape time
				REAL diff_coef=diff_eff*(1e6)/3/B_total; // [eV T⁻¹]
				// = [C⁻¹ kg m² s-2 / (kg C⁻¹ s⁻¹)] = [m² s⁻¹]	
				// to multiply with energy in MeV
				fore_tau=c_r*(1+c_r)*diff_coef/sqr(V_shock); //[s] */
			}		
			REAL losses[7];
			losses[0]= Edot_DSA;
			losses[1]= b_syn;
			losses[2]= b_IC;
			losses[3]= b_br;
			losses[4]= b_coul;
			losses[5]= b_pp;
			losses[6]= b_coulP;
			
			int species;
// arrays for value of cell_n, cell_n+1 and lower edge of cell	
			REAL N[(int)(0.5*(N_OMINT_USER-6)+3)];
			REAL Nnew[(int)(0.5*(N_OMINT_USER-6)+3)];
			REAL Nedge[(int)(0.5*(N_OMINT_USER-6)+3)];
	species=1;
		double deltaT=deltaT_global;

		int new_steps = 1;
		REAL condi = deltaT*Edot(0,losses,nablaV,species);

while (E[1]-E[0] < condi){
//cout << "Edot: " << Edot(0,losses,nablaV,species) << endl;
		//cout << E[1]-E[0] << " and " << condi << endl;
				deltaT=deltaT*0.5;
				new_steps*=2;
			condi = deltaT*Edot(0,losses,nablaV,species);
		}

		//if (deltaT != deltaT_global){
			//cout << "dT=" << deltaT_global << endl;
			//cout << new_steps << " steps with dT = " << deltaT << endl;
			//cout << E[1]-E[0] << " and " << condi << endl;
		//}	
for (int steps=1;steps<=new_steps;steps++){
//if (deltaT != deltaT_global){ cout << "step: " << steps << endl;}
		
			// LOOP over particle species
			for (species=1;species<=2;species++){
				int separator=0;
				if(species==2){
		 			separator=0.5*(N_OMINT_USER-6);
				}
				// initialize fields
				for (bin=0;bin<0.5*(N_OMINT_USER-6);bin++){
					//read in particle density of time n
					N[bin] = gdata.om_user[bin+separator](i,j,k);
					Nnew[bin]=0;
					Nedge[bin]=0;
					// correct for divergent values
					if(N[bin]<=0){
						N[bin]=0;
					}
				}
			//determine highest non-empty energy bin
				int maxi; // index of highest E-bin at t_n
				REAL upper; // energy of next E-bin in t_n+1
				bin=0.5*(N_OMINT_USER-6)-1;
				while(N[bin]<1e-50){
					bin--;
					if(bin==0){
						break;
					}
				}
				maxi=bin;	//highest non empty bin
				upper=E[maxi+1]+deltaT*Edot(maxi+1,losses,
								nablaV,species);
						// value at right border 	

		// PERFORM NUMERICAL SCHEME				
	
				//int ban=0;	// ban := energy bin at time n ; 
				bin=0;		// bin := energy bin at time n+1

				REAL M0 = 0; // := left border of ban projected to n+1:
				REAL M1 = 0; // := right border of ban projected to n+1:
				int cor = 0; // := ban such that its projection contains left border of bin
				double part = 0; // := 1st order approx. of integral over part of ban below Eprime
						

				int low = 0; // ban such that its projection contains left border of bin
				int up = 0; // ban such taht its projection contains right border of bin 
				double partlow = 0; //integral over part of ban projected below left border of bin
				double parthigh = 0; //integral over part of ban projected below right border of bin
								
				while (E[bin]<=upper){

				// first bin has to be treated differenty
					if (bin==0){
						M0 = E[0]+deltaT*Edot(0,losses,nablaV,species);
 						M1 = E[1]+deltaT*Edot(1,losses,nablaV,species);
					} 
				// save identity and integral of previous bin
					else {
						low=cor;
						partlow=part;
					}
				// find ban such that its projection to n+1 contains the left border of bin 
					for(int ban=low;ban<=maxi+1;ban++){
						if ((E[bin+1] >= M0)&&(E[bin+1] < M1)){
							cor=ban;
							break;
						} else {
							M0=M1;  							M1=E[ban+2]+deltaT*Edot(ban+2,losses,nablaV,species);
						}
					}
				// beta := fraction of projected ban below left border of bin
					REAL beta =  (E[bin+1]-M0) / ( M1-M0);
				// Eprime := value in ban such that its projection to n+1 equals left border of bin
					REAL Eprime = E[cor]+beta*deltaE[cor]; 
				// approximate slow of N(E) in low	
					REAL deum = (N[cor]-N[cor-1])/(E[cor+1]-E[cor-1])*2;
					REAL deup = (N[cor+1]-N[cor])/(E[cor+2]-E[cor])*2;
					REAL deuz = (N[cor+1]-N[cor-1])/(E[cor+2]-E[cor-1])*3/2;
					REAL slope = vanLeer(deup,deuz,deum);	
				// approximate left border value of low
					Nedge[cor]=N[cor]-slope*deltaE[cor]*0.5;
				// part:= 1st order approx. of integral over part of ban below Eprime
					part = Nedge[cor]*(Eprime-E[cor]) +
						 0.5*slope*sqr(Eprime-E[cor]);
					if (part<0){
						part=0;
					}
					up=cor;
					parthigh=part;

				if (bin > 0){

					REAL that; // := integral over N(E) at n that is projected to bin =
				// whole of ban MINUS part of low that is not projected to bin
					that= N[low]*deltaE[low]-partlow;
					if (that<= 0){
						that=0;
					}
				// PLUS all bans the projection of which is containded WHOLLY within bin
					int cri=low+1;
					while (cri<up){
						that+=N[cri]*deltaE[cri];
						cri++;
					}
				// PLUS part of up that is no projected to bin+1
					that+=parthigh;
				// SPECIAL CASE: if bin is WHOLLY contained in projection of ban
					if (low==up){
						that=parthigh-partlow;
					}
					if (E[bin+1]>upper){Nnew[bin]= 1e-300;}
				// VALUE of N(E) at bin is integral divided by length of bin		
					Nnew[bin]= that/deltaE[bin];


				// Bohm cutoff		
					if(Ecent[bin]>E_Bohm){  
						Nnew[bin]=0;
					}
		// substract N/T
					//if (Edot_DSA>b_coul){
					if (V_shock >= 0 ){
						//REAL tau=fore_tau*Ecent[bin];
//traditional
						//Nnew[bin]-= Nnew[bin]/T_0*deltaT;
						
			/*
		// add N/T to neighboring cell in WCR direction
						if (direc==1){
							gdata.om_user[bin+separator](i+neigh[0],j+neigh[1],k+neigh[2])
								+=Nnew[bin]/T_0*deltaT;
		//gdata.om_user[N_OMINT_USER-2](i+neigh[0],j+neigh[1],k+neigh[2])=Nnew[1]/T_0*deltaT;
						}else if (direc==2){
							gdata.om_user[bin+separator](i-neigh[0],j-neigh[1],k-neigh[2])
								+=Nnew[bin]/T_0*deltaT;
//gdata.om_user[N_OMINT_USER-2](i-neigh[0],j-neigh[1],k-neigh[2])=Nnew[1]/T_0*deltaT;
						}*/
					
						
					}
		
		
					} //  gt 0
					
					bin++;
				} // end main loop
//traditional
Edot_DSA*=E_zero;
				if(species==1){
					REAL Q_el=inj_el*N_H*1e-6;  // [cm⁻³]
//Q_el=1;
					Nedge[0] =Q_el/(Edot_DSA);//-(b_syn+b_IC)
					//*sqr(E_zero) - b_br*E_zero - b_coul);
//if (Edot_DSA!=0){
//
//}
//else {out2=0;}
				}else if (species==2){
					REAL beta_N0 = sqrt(E_zero*(E_zero+2*M_p))
						/(E_zero+M_p);		
					REAL Q_pr=inj_pr*N_H*1e-6; // [cm⁻³]
//Q_pr=1;
					Nedge[0]=Q_pr/(Edot_DSA);//*Zet-b_pp*E_zero - b_coulP*sqr(beta_N0)/(cube(x_m)+cube(beta_N0)));
				}
				if ((Nedge[0]<0)||(Nedge[0]>1e60)) {
					Nedge[0]=0;
				}
				Nnew[0]=Nedge[0]+(Nedge[1]-Nedge[0])*0.5;
					
				for (bin=0;bin<0.5*(N_OMINT_USER-6);bin++){
					if((Nnew[bin]>0)&&(Nnew[bin]<1e300)) { 
						gdata.om_user[bin+separator](i,j,k)=Nnew[bin];	
						//if (Nnew[bin]>1e-12 && Nnew[bin+1]<1e-12){out2=Ecent[bin];
							//out1=bin;
							//cout << "1: " <<  out2 << endl;
						//}
								
					}
					//if (nablaV > 0){gdata.om_user[bin+separator](i,j,k)=0;}
					//cout << "2: " <<  out1 << endl;
				}
				//cout << "3: " <<  out2 << endl;
				if (species==1){
					gdata.om_user[N_OMINT_USER-4](i,j,k)= B_total;
					gdata.om_user[N_OMINT_USER-1](i,j,k)= 1e-4*std::max(B1,B2);
				}
					
					
			}//species
		} //sub_steps
		gdata.om_user[N_OMINT_USER-5](i,j,k) = V_shock;
		gdata.om_user[N_OMINT_USER-2](i,j,k) = c_r;
	//	gdata.om_user[N_OMINT_USER-1](i,j,k) = out2;	
	
/*if((norm(gdata.getCen_x(i))<gdata.dx[0])&&
(norm(gdata.getCen_y(j))<gdata.dx[0])&&
(norm(gdata.getCen_z(k))<gdata.dx[0])
&& gdata.rank==15)
{cout << gdata.rank << endl;
cout << "electrons: " << gdata.om_user[1](i,j,k) << endl;
cout << "protons: " << gdata.om_user[(int)(0.5*N_OMINT_USER+1)](i,j,k) << endl << endl;
}*/

		}//turn off when in star
		else if( r_sph1 < R_Star1 || r_sph2 < R_Star2){
			for (int bin=0; bin<N_OMINT_USER; bin++){
				gdata.om_user[bin](i,j,k) = -1;	
			}
		}
		else {
			for (int bin=0; bin<N_OMINT_USER; bin++){
				gdata.om_user[bin](i,j,k) = 0;	
			}
		}


	}//i
	}//j
	}//k

//call boundaries 
	for(int bin=0; bin<N_OMINT_USER-6; ++bin) {
		gfunc.boundary(gdata, Problem, gdata.om_user[bin],B,bin);
		
	}
	// LOOP OVER ALL CELLS		
	for(int k=0; k<=gdata.mx[2]; ++k) {
		for(int j=0; j<=gdata.mx[1]; ++j) {
			for(int i=0; i<=gdata.mx[0]; ++i) {
// get distance to stars
				double V_shock = gdata.om_user[N_OMINT_USER-5](i,j,k);
		REAL r_sph1 = sqrt(sqr(gdata.getCen_x(i) - pos_comp1[0]) +
			 	   sqr(gdata.getCen_y(j) - pos_comp1[1]) +
		                   sqr(gdata.getCen_z(k) - pos_comp1[2]));
		REAL r_sph2 = sqrt(sqr(gdata.getCen_x(i) - pos_comp2[0]) +
     	          	           sqr(gdata.getCen_y(j) - pos_comp2[1]) +
     		       	           sqr(gdata.getCen_z(k) - pos_comp2[2]));
				/*if( norm(r_sph1-454.375)<0.1 && norm(r_sph2-985.625)<0.1 && gdata.rank==21){
					cout << V_shock << endl;
					cout << "N: " << gdata.om_user[0](i,j,k) << endl;
					cout << "N: " << gdata.om_user[1](i,j,k) << endl;
					cout << "N: " << gdata.om_user[2](i,j,k) << endl;
				}*/
				double divSix = 1./6;
				
				double c_r = gdata.om_user[N_OMINT_USER-2](i,j,k);
				//double Edot_DSA = gdata.om_user[N_OMINT_USER-4](i,j,k);
				if (V_shock > 0){
				//if ( Edot_DSA > 3.4*gdata.om[q_rho](i,j,k)){
				//diffusion coefficient
					REAL kappa_a = value((char*)"DiffusionConstant"); // [m² s⁻¹]
				// diffusive escape timescale
					double T_0 = c_r*kappa_a/sqr(V_shock); // [s]
					for(int bin=0; bin<N_OMINT_USER-6; ++bin) {
					// add fraction to neighbors
						double present_to_neighbor = divSix*gdata.om_user[bin](i,j,k)/T_0*deltaT_global;
						/*gdata.om_user[bin](i+1,j,k)+=present_to_neighbor;
						gdata.om_user[bin](i-1,j,k)+=present_to_neighbor;
						gdata.om_user[bin](i,j+1,k)+=present_to_neighbor;
						gdata.om_user[bin](i,j-1,k)+=present_to_neighbor;
						gdata.om_user[bin](i,j,k+1)+=present_to_neighbor;
						gdata.om_user[bin](i,j,k-1)+=present_to_neighbor;*/
						gdata.om_user[bin](i,j,k) -= 6*present_to_neighbor;
//if( norm(r_sph1-454.375)<0.1 && norm(r_sph2-985.625)<0.1 && gdata.rank==21){
				/*if (bin<5){
					cout << bin << ": " << gdata.om_user[bin](i,j,k) << " loses " << present_to_neighbor*6 << endl;
					}
}*/
					/*	int direc=1;
						int neigh[3];
						int color;
						neigh[0]=0;
						neigh[1]=0;
						neigh[2]=0;
				// compute gradient of density field
						REAL dTemp_dx=(gdata.om[q_Eges](i+2,j,k)-gdata.om[q_Eges](i-2,j,k))
							/(gdata.getCen_x(i+2)-gdata.getCen_x(i-2));
						REAL dTemp_dy= (gdata.om[q_Eges](i,j+2,k)-gdata.om[q_Eges](i,j-2,k))
							/(gdata.getCen_y(j+2)-gdata.getCen_y(j-2));
						REAL dTemp_dz= (gdata.om[q_Eges](i,j,k+2)-gdata.om[q_Eges](i,j,k-2))
							/(gdata.getCen_z(k+2)-gdata.getCen_z(k-2));
				//and its norm
						REAL dTempnorm = sqrt(sqr(dTemp_dx)+sqr(dTemp_dy)+sqr(dTemp_dz));
			
						//identify neighbor cells along density gradient
						REAL Tempx = norm(dTemp_dx);
						REAL Tempy = norm(dTemp_dy);
						REAL Tempz = norm(dTemp_dz);*/

						/*REAL Tempx = dTemp_dx/dTempnorm;
						REAL Tempy = dTemp_dy/dTempnorm;
						REAL Tempz = dTemp_dz/dTempnorm;

						if (Tempx > 0) { 
							gdata.om_user[bin](i+1,j,k)+=present_to_neighbor*sqr(Tempx);
							gdata.om_user[N_OMINT_USER-1](i+1,j,k)=4;
						}
						else {
							gdata.om_user[bin](i-1,j,k)+=present_to_neighbor*sqr(Tempx);
							gdata.om_user[N_OMINT_USER-1](i-1,j,k)=-4;
						}
						if (Tempy > 0) { 
							gdata.om_user[bin](i,j+1,k)+=present_to_neighbor*sqr(Tempy);
							gdata.om_user[N_OMINT_USER-1](i,j+1,k)=6;
						}
						else {
							gdata.om_user[bin](i,j-1,k)+=present_to_neighbor*sqr(Tempy);
							gdata.om_user[N_OMINT_USER-1](i,j-1,k)=-6;
						}
						if (Tempz > 0) { 
							gdata.om_user[bin](i,j,k+1)+=present_to_neighbor*sqr(Tempz);
							gdata.om_user[N_OMINT_USER-1](i,j,k+1)=8;
						}
						else {
							gdata.om_user[bin](i,j,k-1)+=present_to_neighbor*sqr(Tempz);
							gdata.om_user[N_OMINT_USER-1](i,j,k-1)=-8;
						}*/
		

					/*	if (Tempx>=Tempy && Tempx>=Tempz){ // neighbors along x
							neigh[0]=1; 
							color = 4;
						}else if (Tempy>=Tempx && Tempy>=Tempz){
							neigh[1]=1; 
							color = 6;
						}else if (Tempz>=Tempx && Tempz>=Tempy){
							neigh[2]=1; 
							color = 8;
						}
						// get density of neighbor cells
						int thick=2;
						REAL Temp_neigh1 = 
						gdata.om[q_Eges](i+thick*neigh[0],j+thick*neigh[1],k+thick*neigh[2]);
						REAL Temp_neigh2 = 
						gdata.om[q_Eges](i-thick*neigh[0],j-thick*neigh[1],k-thick*neigh[2]);
						if (Temp_neigh2 > Temp_neigh1){direc=2;}
						if (direc==1){
							gdata.om_user[bin](i+neigh[0],j+neigh[1],k+neigh[2])
								+=6*present_to_neighbor;
						//=gdata.om_user[bin](i,j,k);
							gdata.om_user[N_OMINT_USER-1](i+neigh[0],j+neigh[1],k+neigh[2])=color;
		
						}else if (direc==2){
							gdata.om_user[bin](i-neigh[0],j-neigh[1],k-neigh[2])
								+=6*present_to_neighbor;
						//=gdata.om_user[bin](i,j,k);
							gdata.om_user[N_OMINT_USER-1](i-neigh[0],j-neigh[1],k-neigh[2])=color;

						}*/
					
				
					
					}
				}
			}
		}
	}
//call boundaries again
	for(int bin=0; bin<N_OMINT_USER-6; ++bin) {
		gfunc.boundary(gdata, Problem, gdata.om_user[bin],B,bin);
	}
	

#endif
	return 0;	
	
}


REAL semiLagra::norm(REAL x){
	return sqrt(x*x);
}

REAL semiLagra::vanLeer(REAL &up, REAL &uc, REAL &um){
	return (std::max((up*um),0.))/(uc+1.e-18);
}

REAL semiLagra::DiLog(REAL x){
	REAL preis = 0;	
	if (x>-1){
		for (int i = 1; i<=100; i++){
			preis += pow(x,i)/sqr(i);
		}
	}else if (x<-1){
		for (int i = 1; i<=100; i++){
			preis += pow(1/(1-x),i)/sqr(i);
		}
		preis = -sqr(M_PI)/6 - 0.5*log(1-x)*( 2*log(-x)-log(1-x) ) + preis;
	}else cout << "Error: x= " << x << endl;
	
	return preis;
}

//determination of MAGNETIC FIELD STRENGTH following Usov and Melrose 1992
REAL semiLagra::findB(int which, REAL dist){ // called for every grid point and timestep

	REAL Bfield=0;	// [G]
	if (which==1){
		if (dist < rA_Star1) {	
			//cout << "OB wind suppressed due to too strong WR wind !" << endl;
			Bfield = B_surface_Star1*cube(R_Star1/dist);
		}else if ( (dist >= rA_Star1) && (dist< x_toro_Star1) ){
			// cout << "mainly radial field, shock parallel" << endl;
			Bfield = B_surface_Star1*cube(R_Star1)/rA_Star1/sqr(dist);
      		}else if ( (dist >= rA_Star1) && (dist >= x_toro_Star1) ){
			 //cout << "mainly toroidal magnetic field, perpendicular shock" << endl;
			Bfield = B_surface_Star1*cube(R_Star1)/rA_Star1/dist/x_toro_Star1; 
			Bfield = B_surface_Star1*cube(R_Star1)/rA_Star1/sqr(dist);
		}
	}else if (which==2){
		if (dist < rA_Star2) {	
			//cout << "OB wind suppressed due to too strong WR wind !" << endl;
			Bfield = B_surface_Star2*cube(R_Star2/dist);
		}else if ( (dist >= rA_Star2) && (dist< x_toro_Star2) ){
			 //cout << "mainly radial field, shock parallel" << endl;
			Bfield = B_surface_Star2*cube(R_Star2)/rA_Star2/sqr(dist);
      		}else if ( (dist >= rA_Star2) && (dist >= x_toro_Star2) ){
			// cout << "mainly toroidal magnetic field, perpendicular shock" << endl;
			Bfield = B_surface_Star2*cube(R_Star2)/rA_Star2/dist/x_toro_Star2; 
			Bfield = B_surface_Star2*cube(R_Star2)/rA_Star2/sqr(dist);
		}
	}		
	return Bfield; // [G]
}

REAL semiLagra::Edot(int james,REAL *losses, REAL nablaV, int species){
	
	REAL E = Ecent[james];
	REAL KN = full_KN[james];
	REAL Edot1 = 0;
//traditional
	if(species==1){
		Edot1= losses[0]*E-(losses[1]+losses[2]*KN)*sqr(E)-losses[3]*E-losses[4];
	}
	else if (species==2){
		double beta_N = sqrt(E*(E+2*M_p))/(E+M_p);
		Edot1 = Zet*losses[0]*E-losses[5]*E
			 - losses[6]*sqr(beta_N)/(cube(x_m)+cube(beta_N));
	}
//NEW!!!!!!!!!!!!!!!!!!!!!!
	if(nablaV > 0){
		Edot1 -= 0.3333*E* nablaV *c_sound / RSun; // [MeV s⁻¹]
	}			 

	return Edot1;
}

REAL semiLagra::EccentricAnomaly(REAL normtime, REAL psi)
{
/*
	solve omega*t = psi - e*sin(psi);
	using the Newton-Raphson method, i.e. we solve f(u)=0
	with the formula unew=uold-f(uold)/f'(uold)
	We have 0=f(u)= psi - e*sin(psi) - omega*t
	f'(u) = 1 - e*cos(psi)
	*/

	REAL u= psi;
	REAL time = normtime*RSun/c_sound;
	REAL omega1 = sqrt(G_Grav*MSun*(M_Star1 + M_Star2)/cube(semimajor*RSun));
	int manfred;
	if (apastron==1){manfred=M_PI/omega1;}
	else {manfred=0;}

/* 50 iterations should be enough */

	for (int n=1; n<50; n++) {
		REAL unew = u - ( (u - eccentricity*sin(u) - omega1*(time+manfred))/(1 - eccentricity*cos(u)) );
		// cout << " fu: " << psi << " " << psinew << " ";
		//cout << endl;
		u = unew;
	}
	
	// Multiply by normalised cs2 to find result.
	return u;

}

int semiLagra::getEnergy(){
	int bin;
	REAL deltalog;
						
	for (bin=0;bin<=0.5*N_OMINT_USER-2;bin++){
		
		//non-linear grid
		deltalog=(log10(Emax)-log10(Emin))/(0.5*(REAL)N_OMINT_USER-3);
		E[bin]= log10(Emin)+bin*deltalog;
		E[bin]= pow(10,E[bin]);
		deltaE[bin]= log10(Emin)+(bin+1)*deltalog;
		deltaE[bin]= pow(10,deltaE[bin]) - E[bin];
			
		Ecent[bin]= 0.5*(log10(Emin)+(bin+1)*deltalog+log10(E[bin]));
		Ecent[bin]= pow(10,Ecent[bin]);
		REAL b= 4*Ecent[bin]*(eps_T/1e6)/sqr(0.511);
		if (b>=1e-3){
			full_KN[bin] =  (0.5*b+6+6/b)*log(1+b);
			full_KN[bin] += -(0.916666666*cube(b)+6*sqr(b)+9*b+4)/sqr(1+b);
			full_KN[bin] += -2 + 2*DiLog(-b);
			full_KN[bin] *= 9/cube(b);
	
		}
		else {full_KN[bin] = 1;}
	
		//full_KN[bin] = 1/pow(1+b,1.5); // TURN ON for Moderski
		//cout << bin << " E: " << E[bin] << endl;
		//cout << " [" << setprecision(2) << E[bin-1] << ", " << setprecision(2) << E[bin] << "] MeV ";
	}
	return 0;
}
void semiLagra::writesomething(Data &gdata){
if(gdata.rank==0){
			cout << "alfven1: " << rA_Star1 << " and alfven2: " << rA_Star2 << endl;
			cout << "x_toro1: " << x_toro_Star1 << "and x_toro2: " << x_toro_Star2 << endl;
		}
return;
}
	

	
	


