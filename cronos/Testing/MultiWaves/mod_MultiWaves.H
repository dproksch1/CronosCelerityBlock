#ifndef CRONOS_MOD_MULTIWAVES_H
#define CRONOS_MOD_MULTIWAVES_H = 1

#include "problem.H"
#include "data.H"
#include "gridgen.H"
#include <string>

class MultiWaves: public ProblemType {
public:
	MultiWaves (const Data &);
	virtual void init_fields (Data &, int [], int []);
	void src_User (Data &gdata,
	               NumMatrix<REAL,3> nom     [N_OMINT     ],
	               NumMatrix<REAL,3> nom_user[N_OMINT_USER]);
private:
//  It all happens in the xz plane! (y direction invariant)
//   k and Bo along x, dB and dv along z
	
	REAL B0mag;      // background field magnitude
	REAL B1mag;      // amplitude of magnetic disturbance
	REAL dens0;      // base density
	REAL kx;           // wave number
	REAL om_re, om_im; // Re(omega), Im(omega)
	bool isothermal;   // (maybe) for later use...
};


MultiWaves::MultiWaves (const Data &gdata) : ProblemType(gdata) {

#if (FLUID_TYPE != CRONOS_MULTIFLUID)
	cerr << endl;
	cerr << "!---------------------------------------------!"<< endl;
	cerr << "! ERROR: Need to use FLUID_TYPE CRONOS_MULTIFLUID    !"<< endl;
	cerr << "!    --> correct in constants.H! Exiting...   !"<< endl;
	cerr << "!---------------------------------------------!"<< endl;
	cerr << endl;
	exit(2);
#endif

	if(ENERGETICS == FULL) {
		name = " Multifluid waves test";
		isothermal = false;
	} else {
		name = " Isothermal multifluid waves test";
		if(gamma > 1.) {
			cerr << " Has to be isothermal ";
			exit(2);
		}
		isothermal = true;
	}
	
#if (ENERGETICS == FULL)
	cerr << "Only isothermal, sorry!" << endl;
	exit(1);
#endif
	REAL xlen = value((char*)"xe")-value((char*)"xb"); // w/o pi factor
	kx = int(2./xlen * 100.)/100.;  // wavenumber (rounded to 0.01)

	B0mag = value((char*)"B0mag");
	B1mag = value((char*)"B1mag");
	dens0 = value((char*)"dens0");
	
	ifstream omega_table ("tab_k-omRE-omIM.dat");
	if (omega_table.fail()) {
		cerr << "ERROR: cannot access read-in table. Aborting." << endl;
		exit(1);
	}
	
	REAL a = 0.;  // default
	REAL t_a;     // temp. read-in var
	std::istringstream lin;
	for (std::string line; std::getline(omega_table, line); ) {
		
		std::size_t beg_comment = line.find_first_of ("#");
		if (beg_comment != std::string::npos)   // ignore comments
			line.resize (beg_comment);           // (starting with '#')
		
		if (beg_comment > 0) {
			lin.clear();
			lin.str (line);
			if (lin >> t_a >> om_re >> om_im) {
				a = t_a;
				if (t_a == kx) break; // look-up value found.
			}
			else
				cerr << "WARNING: failed to decode line '" << line << "'" << endl;
		}
		if (a == 100.) {   // end of file reached
			cerr << "ERROR: wavenumber a=" << kx << " not found in table." << endl;
			exit(1);
		}
	}
}


void MultiWaves::init_fields (Data &gdata, int ibeg[3], int iend[3])
{	
	REAL xi_n = 0.5;       // neutral fraction of density (hard-coded)
	int numFluids = gdata.fluids->get_numFluids();
	q_Bx = gdata.fluids->get_q_Bx();
	q_By = gdata.fluids->get_q_By();
	q_Bz = gdata.fluids->get_q_Bz();
	
	for (int iFluid = 0; iFluid < numFluids; ++iFluid) {
		int fluidType = gdata.fluids->get_fluidType(iFluid);
		cout << " Init for fluid " << iFluid << " of type " << fluidType << endl;

		q_rho = gdata.fluids->fluids[iFluid].get_q_rho_global();
		q_sz  = gdata.fluids->fluids[iFluid].get_q_sz_global();
		
		// coefficents for cos,sin
		REAL coco[] = {om_re/kx,               // ionized
		 					om_re/kx*(om_im+1.) };  // neutral
		REAL cosi[] = {om_im/kx,               // ionized
							(om_re*om_re-om_im*om_im)/(2.*kx) -om_im/kx -kx};
		// REAL coco[] = {1, 1};  // hack
		// REAL cosi[] = {0, 0};
		
		for (int k = ibeg[2]; k <= iend[2]; k++){
			//REAL zz = gdata.getCen_z(k);
			for (int j = ibeg[1]; j <= iend[1]; j++){
				//REAL yy = gdata.getCen_y(j);
				for (int i = ibeg[0]; i <= iend[0]; i++){
					REAL xx = gdata.get_x (i, 0);
					REAL xr = gdata.get_x (i,+1);
					
					if (iFluid == 0) {  // ionized
						gdata.om[q_rho](i,j,k) = (1.-xi_n)*dens0;
						gdata.om[q_sz ](i,j,k) = -B1mag*(coco[0]*cos(kx*xx)-
																	cosi[0]*sin(kx*xx));
					} else {            // neutral
						gdata.om[q_rho](i,j,k) =     xi_n *dens0;
						gdata.om[q_sz ](i,j,k) = -B1mag*(coco[1]*cos(kx*xx)+
																	cosi[1]*sin(kx*xx));
					}
					if (fluidType == CRONOS_MHD) {
						gdata.om[q_Bx](i,j,k) = B0mag;
						gdata.om[q_Bz](i,j,k) = B1mag*cos(kx*xr);
					}
				}
			}
		}
	}
}


void MultiWaves::src_User (Data &gdata,
									NumMatrix<REAL,3> nom     [N_OMINT     ],
									NumMatrix<REAL,3> nom_user[N_OMINT_USER])
{
	int q_sz_i = gdata.fluids->fluids[0].get_q_sz_global();  // Vz_ion
	int q_sz_n = gdata.fluids->fluids[1].get_q_sz_global();  // Vz_ntr
	
	for (int k = 0; k <= gdata.mx[2]; k++)
		for (int j = 0; j <= gdata.mx[1]; j++)
			for (int i = 0; i <= gdata.mx[0]; i++) {
				REAL ui_un = gdata.om[q_sz_i](i,j,k) - gdata.om[q_sz_n](i,j,k);
				nom[q_sz_i](i,j,k) -= -ui_un;
				nom[q_sz_n](i,j,k) -= +ui_un;
			}
}

#endif // [CRONOS_MOD_MULTIWAVES_H]
