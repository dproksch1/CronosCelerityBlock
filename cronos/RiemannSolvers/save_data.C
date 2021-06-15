#include "save_data.H"
#include <iostream>

using namespace std;

Saves::Saves(const Data& gdata) {
#if (FLUID_TYPE != CRONOS_HYDRO)

#if (FLUID_TYPE == CRONOS_MHD)
	q_rho = gdata.fluid.get_q_rho();
	q_sx = gdata.fluid.get_q_sx();
	q_sy = gdata.fluid.get_q_sy();
	q_sz = gdata.fluid.get_q_sz();
	q_Bx = gdata.fluid.get_q_Bx();
	q_By = gdata.fluid.get_q_By();
	q_Bz = gdata.fluid.get_q_Bz();
#elif (FLUID_TYPE == CRONOS_MULTIFLUID)
	int i_magFluid = gdata.fluids->get_i_magFluid();
	q_rho = gdata.fluids->fluids[i_magFluid].get_q_rho();
	q_sx = gdata.fluids->fluids[i_magFluid].get_q_sx();
	q_sy = gdata.fluids->fluids[i_magFluid].get_q_sy();
	q_sz = gdata.fluids->fluids[i_magFluid].get_q_sz();
	q_Bx = gdata.fluids->get_q_Bx();
	q_By = gdata.fluids->get_q_By();
	q_Bz = gdata.fluids->get_q_Bz();
#endif

	int num(0);
#if CT_TYPE == CONSISTENT
	num = 7;
#elif CT_TYPE == STONE
	num = 2;
#if STONE_TYPE == STONE_CENTRE
	// We need to save the mass flux additionally
	num = 3;
#endif
#endif

	save_field_x = new NumMatrix<REAL,3> [num];
	save_field_y = new NumMatrix<REAL,3> [num];
	save_field_z = new NumMatrix<REAL,3> [num];

	for (int q = 0; q < num; ++q) {
		save_field_x[q].resize(Index::set(-2,-2,-2),
		                       Index::set(gdata.mx[0]+1,gdata.mx[1]+1,gdata.mx[2]+1));
		save_field_y[q].resize(Index::set(-2,-2,-2),
		                       Index::set(gdata.mx[0]+1,gdata.mx[1]+1,gdata.mx[2]+1));
		save_field_z[q].resize(Index::set(-2,-2,-2),
		                       Index::set(gdata.mx[0]+1,gdata.mx[1]+1,gdata.mx[2]+1));

		// Initialisation of values
		save_field_x[q].set_constVal(0.);
		save_field_y[q].set_constVal(0.);
		save_field_z[q].set_constVal(0.);
		
	}

#endif
}

Saves::~Saves() {
#if (FLUID_TYPE != CRONOS_HYDRO)
	delete [] save_field_x;
	delete [] save_field_y;
	delete [] save_field_z;
#endif
}

void Saves::save_vars(const Data &gdata, fields_1D &fields, const int &dir,
                      const int &j, const int &k) {
#if (FLUID_TYPE != CRONOS_HYDRO)

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	int i_magFluid = gdata.fluids->get_i_magFluid();
	q_rho = gdata.fluids->get_q_rho(i_magFluid);
	q_sx = gdata.fluids->get_q_sx(i_magFluid);
	q_sy = gdata.fluids->get_q_sy(i_magFluid);
	q_sz = gdata.fluids->get_q_sz(i_magFluid);
	q_Bx = gdata.fluids->get_q_Bx(i_magFluid);
	q_By = gdata.fluids->get_q_By(i_magFluid);
	q_Bz = gdata.fluids->get_q_Bz(i_magFluid);
#endif


	if(dir == 0) {
		for (int i = -2; i <= gdata.mx[dir]+1; ++i){

#if CT_TYPE == CONSISTENT

#ifndef CRONOS_SAVEMEM
			// Saving momentum reconstruction
			save_field_x[0](i,j,k) = fields.deriv[q_sx](i);
			save_field_x[1](i,j,k) = fields.deriv[q_sy](i);
			save_field_x[2](i,j,k) = fields.deriv[q_sz](i);

			// Saving magnetic induction
			save_field_x[3](i,j,k) = fields.deriv[q_By](i);
			save_field_x[4](i,j,k) = fields.deriv[q_Bz](i);
#endif
			// Saving characteristic velocities
			save_field_x[5](i,j,k) = fields.v_ch_p(i);
			save_field_x[6](i,j,k) = fields.v_ch_m(i);

#elif CT_TYPE == STONE

			save_field_x[0](i,j,k) = fields.flux[q_By](i);
			save_field_x[1](i,j,k) = fields.flux[q_Bz](i);

#if STONE_TYPE == STONE_CENTRE

			save_field_x[2](i,j,k) = fields.flux[q_rho](i);

#endif

#endif

		}
	} else if (dir == 1) {
		for (int i = -2; i <= gdata.mx[dir]+1; ++i){

#if CT_TYPE == CONSISTENT

#ifndef CRONOS_SAVEMEM
			// Saving momentum reconstruction
			save_field_y[0](j,i,k) = fields.deriv[q_sx](i);
			save_field_y[1](j,i,k) = fields.deriv[q_sy](i);
			save_field_y[2](j,i,k) = fields.deriv[q_sz](i);

			// Saving magnetic induction
			save_field_y[3](j,i,k) = fields.deriv[q_Bx](i);
			save_field_y[4](j,i,k) = fields.deriv[q_Bz](i);
#endif
			// Saving characteristic velocities
			save_field_y[5](j,i,k) = fields.v_ch_p(i);
			save_field_y[6](j,i,k) = fields.v_ch_m(i);

#elif CT_TYPE == STONE

			save_field_y[0](j,i,k) = fields.flux[q_Bx](i);
			save_field_y[1](j,i,k) = fields.flux[q_Bz](i);

#if STONE_TYPE == STONE_CENTRE

			save_field_y[2](j,i,k) = fields.flux[q_rho](i);

#endif

#endif

		}
	} else if (dir == 2) {
		for (int i = -2; i <= gdata.mx[dir]+1; ++i){

#if CT_TYPE == CONSISTENT

#ifndef CRONOS_SAVEMEM
			// Saving momentum reconstruction
			save_field_z[0](j,k,i) = fields.deriv[q_sx](i);
			save_field_z[1](j,k,i) = fields.deriv[q_sy](i);
			save_field_z[2](j,k,i) = fields.deriv[q_sz](i);

			// Saving magnetic induction
			save_field_z[3](j,k,i) = fields.deriv[q_Bx](i);
			save_field_z[4](j,k,i) = fields.deriv[q_By](i);
#endif
			// Saving characteristic velocities
			save_field_z[5](j,k,i) = fields.v_ch_p(i);
			save_field_z[6](j,k,i) = fields.v_ch_m(i);

#elif CT_TYPE == STONE

			save_field_z[0](j,k,i) = fields.flux[q_Bx](i);
			save_field_z[1](j,k,i) = fields.flux[q_By](i);

#if STONE_TYPE == STONE_CENTRE

			save_field_z[2](j,k,i) = fields.flux[q_rho](i);

#endif

#endif

		}
	}
#endif
}


void Saves::retrieve(const Data &gdata, fields_2D &fields,
                     const int &dir, const int &layer) {

#if (FLUID_TYPE != CRONOS_HYDRO)

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	int i_magFluid = gdata.fluids->get_i_magFluid();
	q_rho = gdata.fluids->get_q_rho(i_magFluid);
	q_sx = gdata.fluids->get_q_sx(i_magFluid);
	q_sy = gdata.fluids->get_q_sy(i_magFluid);
	q_sz = gdata.fluids->get_q_sz(i_magFluid);
	q_Bx = gdata.fluids->get_q_Bx(i_magFluid);
	q_By = gdata.fluids->get_q_By(i_magFluid);
	q_Bz = gdata.fluids->get_q_Bz(i_magFluid);
#endif

#ifndef CRONOS_SAVEMEM
	if(dir == 0) {
		for (int k = -2; k <= gdata.mx[2]+1; ++k){
			for (int j = -2; j <= gdata.mx[1]+1; ++j){

#if CT_TYPE == CONSISTENT
	
				fields.dvydx[0](j,k) = save_field_y[1](layer,j,k); // dvydy
				fields.dvydx[1](j,k) = save_field_z[1](layer,j,k); // dvydz
	
				fields.dvzdx[0](j,k) = save_field_y[2](layer,j,k); // dvzdy
				fields.dvzdx[1](j,k) = save_field_z[2](layer,j,k); // dvzdz
	
				fields.dBydx(j,k) = save_field_z[4](layer,j,k); // dBydz
				fields.dBzdx(j,k) = save_field_y[4](layer,j,k); // dBzdy
	
				fields.v_ch_p[0](j,k) = save_field_y[5](layer,j,k);
				fields.v_ch_p[1](j,k) = save_field_z[5](layer,j,k);
	
				fields.v_ch_m[0](j,k) = save_field_y[6](layer,j,k);
				fields.v_ch_m[1](j,k) = save_field_z[6](layer,j,k);
	
#elif CT_TYPE == STONE

				fields.fluxSN[q_Bz](j,k) = save_field_y[1](layer,j,k);
				fields.fluxBT[q_By](j,k) = save_field_z[1](layer,j,k);

#if STONE_TYPE == STONE_CENTRE
				fields.fluxSN[q_rho](j,k) = save_field_y[2](layer,j,k);
				fields.fluxBT[q_rho](j,k) = save_field_z[2](layer,j,k);
#endif

#endif

			}
		}
	} else if (dir == 1) {
		for (int k = -2; k <= gdata.mx[2]+1; ++k){
			for (int i = -2; i <= gdata.mx[0]+1; ++i){
	
#if CT_TYPE == CONSISTENT

				fields.dvxdx[0](k,i) = save_field_z[0](i,layer,k); // dvxdz
				fields.dvxdx[1](k,i) = save_field_x[0](i,layer,k); // dvxdx
	
				fields.dvzdx[0](k,i) = save_field_z[2](i,layer,k); // dvzdz
				fields.dvzdx[1](k,i) = save_field_x[2](i,layer,k); // dvzdx
	
				fields.dBxdx(k,i) = save_field_z[3](i,layer,k); // dBxdz
				fields.dBzdx(k,i) = save_field_x[4](i,layer,k); // dBzdx
	
				fields.v_ch_p[0](k,i) = save_field_z[5](i,layer,k);
				fields.v_ch_p[1](k,i) = save_field_x[5](i,layer,k);
	
				fields.v_ch_m[0](k,i) = save_field_z[6](i,layer,k);
				fields.v_ch_m[1](k,i) = save_field_x[6](i,layer,k);
	
#elif CT_TYPE == STONE

				fields.fluxWE[q_Bz](k,i) = save_field_x[1](i,layer,k);
				fields.fluxBT[q_Bx](k,i) = save_field_z[0](i,layer,k);

#if STONE_TYPE == STONE_CENTRE
				fields.fluxWE[q_rho](k,i) = save_field_x[2](i,layer,k);
				fields.fluxBT[q_rho](k,i) = save_field_z[2](i,layer,k);
#endif


#endif

			}
		}
	} else {
		for (int j = -2; j <= gdata.mx[1]+1; ++j){
			for (int i = -2; i <= gdata.mx[0]+1; ++i){
	
#if CT_TYPE == CONSISTENT

				fields.dvxdx[0](i,j) = save_field_x[0](i,j,layer); // dvxdx
				fields.dvxdx[1](i,j) = save_field_y[0](i,j,layer); // dvxdy
	
				fields.dvydx[0](i,j) = save_field_x[1](i,j,layer); // dvydx
				fields.dvydx[1](i,j) = save_field_y[1](i,j,layer); // dvydy
	
				fields.dBxdx(i,j) = save_field_y[3](i,j,layer); // dBxdy
				fields.dBydx(i,j) = save_field_x[3](i,j,layer); // dBydx
	
				fields.v_ch_p[0](i,j) = save_field_x[5](i,j,layer);
				fields.v_ch_p[1](i,j) = save_field_y[5](i,j,layer); 
	
				fields.v_ch_m[0](i,j) = save_field_x[6](i,j,layer);
				fields.v_ch_m[1](i,j) = save_field_y[6](i,j,layer); 
	
#elif CT_TYPE == STONE

				fields.fluxWE[q_By](i,j) = save_field_x[0](i,j,layer);
				fields.fluxSN[q_Bx](i,j) = save_field_y[0](i,j,layer);

#if STONE_TYPE == STONE_CENTRE
				fields.fluxWE[q_rho](i,j) = save_field_x[2](i,j,layer);
				fields.fluxSN[q_rho](i,j) = save_field_y[2](i,j,layer);
#endif

#endif

			}
		}
    
	}
#endif
#endif
}

