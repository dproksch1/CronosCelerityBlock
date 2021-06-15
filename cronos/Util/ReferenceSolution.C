#include "ReferenceSolution.H"
#include "Hdf5File_cbase.H"
#include "CException.H"
#include <iostream>
#include <sys/stat.h>
#ifdef parallel
#include "mpi.h"
#endif

using namespace std;

ReferenceSolution::ReferenceSolution(string full_path, int rank) {
	this->full_path = full_path;
	this->mpi_rank = rank;
	this->mx_ref.resize(3);

	this->i_refNearestL.resize(3);
	this->dist_refPoints.resize(3);
	this->dist_fromRefL.resize(3);
	this->dist_fromRefR.resize(3);

	this->do_tilt = false;
	this->theta_tilt = 0.;

	load_file();
}


void ReferenceSolution::set_tiltAngle(double tilt_angle) {
	this->theta_tilt = tilt_angle;
	this->do_tilt = true;
}

#ifdef parallel
void ReferenceSolution::broadcast_map_entry(std::string &name_property, float &value_property) {
	int name_size = name_property.size();
	MPI_Bcast(&name_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(this->mpi_rank != 0) {
		name_property.resize(name_size);
	}
	MPI_Bcast(const_cast<char*>(name_property.data()), name_size, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast(&value_property, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
}



void ReferenceSolution::broadcast_data() {

	// Start by broadcasting the properties-map from rank 0:
	int map_size = map_reference_properties.size();
	MPI_Bcast(&map_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(this->mpi_rank==0) {
		for (auto &map_entry: map_reference_properties) {

			string name_entry = map_entry.first;
			float val_entry = map_entry.second;

			broadcast_map_entry( name_entry, val_entry );
		}
	} else {
		for(int i_entry=0; i_entry<map_size; ++i_entry) {

			string name_entry;
			float val_entry;

			broadcast_map_entry( name_entry, val_entry );

			map_reference_properties[ name_entry ] = val_entry;

		}
	}

	// Synchronise information on the grid
	MPI_Bcast(&grid_type_rad, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(grid_type_rad == 6) {
		MPI_Bcast(&grid_beg_rad, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&grid_inc_rad, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&grid_del_rad, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}

	MPI_Bcast(&grid_type_theta, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(grid_type_theta == 0) {
		MPI_Bcast(&grid_beg_theta, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&grid_end_theta, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&grid_del_theta, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}

	// Next, we synchronise the fluid data
	MPI_Bcast(mx_ref, 3, MPI_INT, 0, MPI_COMM_WORLD);


	// First, we need to resize all arrays accordingly:
	if(mpi_rank != 0) {
		radCen.resize(mx_ref[0]+1);
		radL.resize(mx_ref[0]+2);
		thetaCen_ref.resize(mx_ref[1]+1);
		thetaL_ref.resize(mx_ref[1]+2);

		rho_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0], mx_ref[1], mx_ref[2]));
		vr_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0], mx_ref[1], mx_ref[2]));
		vtheta_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0], mx_ref[1], mx_ref[2]));
		vphi_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0], mx_ref[1], mx_ref[2]));
		Ar_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0], mx_ref[1]+1, mx_ref[2]+1));
		Atheta_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0]+1, mx_ref[1], mx_ref[2]+1));
		Aphi_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0]+1, mx_ref[1]+1, mx_ref[2]));
		dArdt_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0], mx_ref[1]+1, mx_ref[2]+1));
		dAthetadt_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0]+1, mx_ref[1], mx_ref[2]+1));
		dAphidt_ref.resize(Index::set(0,0,0),
				Index::set(mx_ref[0]+1, mx_ref[1]+1, mx_ref[2]));

	}

	MPI_Barrier(MPI_COMM_WORLD);


	// Then, we transfer the array data
	MPI_Bcast(radCen, mx_ref[0]+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(radL, mx_ref[0]+2, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(thetaCen_ref, mx_ref[1]+1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(thetaL_ref, mx_ref[1]+2, MPI_FLOAT, 0, MPI_COMM_WORLD);

	int size = (mx_ref[0]+1)*(mx_ref[1]+1)*(mx_ref[2]+1);
	MPI_Bcast(rho_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vr_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vtheta_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(vphi_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);

	size = (mx_ref[0]+1)*(mx_ref[1]+2)*(mx_ref[2]+2);
	MPI_Bcast(Ar_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(dArdt_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);

	size = (mx_ref[0]+2)*(mx_ref[1]+1)*(mx_ref[2]+2);
	MPI_Bcast(Atheta_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(dAthetadt_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);

	size = (mx_ref[0]+2)*(mx_ref[1]+2)*(mx_ref[2]+1);
	MPI_Bcast(Aphi_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(dAphidt_ref, size, MPI_FLOAT, 0, MPI_COMM_WORLD);


}
#endif


bool ReferenceSolution::file_exists(const std::string& filename) {
  struct stat buffer;
  return (stat (filename.c_str(), &buffer) == 0);
}


int ReferenceSolution::load_file() {

	if( this->mpi_rank==0 ) {

		if(!file_exists(this->full_path)) {
			string message = "Reference file " + this->full_path + " does not exist ";
			throw CException("message");
		}

		cout << " Trying to read from file " << full_path << endl;
		// Do not open the file in parallel mode
		Hdf5iStream h5in(this->full_path, mpi_rank, false);

		// Load position arrays
		h5in.ReadNumArray("r_Cen", radCen);
		h5in.ReadNumArray("r_L", radL);
		h5in.ReadNumArray("theta_Cen", thetaCen_ref);
		h5in.ReadNumArray("theta_L", thetaL_ref);

		// Now load 3D data
		h5in.Read3DMatrix("rho", rho_ref, true);
		h5in.Read3DMatrix("v_x", vr_ref, true);
		h5in.Read3DMatrix("v_y", vtheta_ref, true);
		h5in.Read3DMatrix("v_z", vphi_ref, true);
		h5in.Read3DMatrix("A_x", Ar_ref, true);
		h5in.Read3DMatrix("A_y", Atheta_ref, true);
		h5in.Read3DMatrix("A_z", Aphi_ref, true);
		h5in.Read3DMatrix("dAxdt", dArdt_ref, true);
		h5in.Read3DMatrix("dAydt", dAthetadt_ref, true);
		h5in.Read3DMatrix("dAzdt", dAphidt_ref, true);

//		cout << " Reading " << endl;
//		for(int ir=0; ir<350; ir++) {
//				cout << ir << " " << vr_ref(ir,10,0) << " " << vr_ref.getHigh(0) << endl;
//			}

		// Determine size of arrays
		for(int dir=0; dir<3; ++dir) {
			mx_ref[dir] = rho_ref.getHigh(dir);
		}

		if(h5in.doesSubGroupExist("GridPars")) {
			string group_name = "Data/GridPars";
			// radial grid
			h5in.ReadGlobalAttr(group_name, "grid_type_x", grid_type_rad);


			// So far only for grid type 6:
			if(grid_type_rad == 6) {
				h5in.ReadGlobalAttr(group_name, "xb_x", grid_beg_rad);
				h5in.ReadGlobalAttr(group_name, "cell_change_x", grid_inc_rad);
				h5in.ReadGlobalAttr(group_name, "ref_size_x", grid_del_rad);
			}

			// theta grid
			h5in.ReadGlobalAttr(group_name, "grid_type_y", grid_type_theta);

			// So far only for grid type 0
			if(grid_type_theta == 0) {
				h5in.ReadGlobalAttr(group_name, "xb_y", grid_beg_theta);
				float gridEndT_Star1;
				h5in.ReadGlobalAttr(group_name, "xe_y", grid_end_theta);
				int Ny_ref;
				h5in.ReadGlobalAttr(group_name, "N_y", Ny_ref);
				grid_del_theta =  (grid_end_theta-grid_beg_theta)/(1.*Ny_ref);
				if(Ny_ref != mx_ref[1]+1) {
					cerr << " Error: discrepancy for for number of grid points in theta direction " << endl;
				}
			}


			// If stellar properties are stores: load corresponding data
			if(h5in.doesSubGroupExist("StellarProperties")) {
				string group_name = "Data/StellarProperties";
				float k_CAK, alpha_CAK, v_terminal_norm, v_terminal_SI;
				float StellarMass, StellarTemperature, MassLossRate_SI, MassLossRate_norm, StellarRadius_SI, StellarRadius_norm;
				h5in.ReadGlobalAttr(group_name, "k_CAK", k_CAK);
				h5in.ReadGlobalAttr(group_name, "alpha_CAK", alpha_CAK);
				h5in.ReadGlobalAttr(group_name, "v_terminal", v_terminal_SI);
				h5in.ReadGlobalAttr(group_name, "v_terminal_norm", v_terminal_norm);
				h5in.ReadGlobalAttr(group_name, "StellarMass", StellarMass);
				h5in.ReadGlobalAttr(group_name, "StellarTemperature", StellarTemperature);
				h5in.ReadGlobalAttr(group_name, "StellarRadius", StellarRadius_SI);
				h5in.ReadGlobalAttr(group_name, "StellarRadius_norm", StellarRadius_norm);
				h5in.ReadGlobalAttr(group_name, "MasslossRate", MassLossRate_SI);
				h5in.ReadGlobalAttr(group_name, "MasslossRate_norm", MassLossRate_norm);

				map_reference_properties["k_CAK"] = k_CAK;
				map_reference_properties["alpha_CAK"] = alpha_CAK;
				map_reference_properties["v_terminal_SI"] = v_terminal_SI;
				map_reference_properties["v_terminal_norm"] = v_terminal_norm;
				map_reference_properties["StellarMass"] = StellarMass;
				map_reference_properties["StellarTemperature"] = StellarTemperature;
				map_reference_properties["StellarRadius_SI"] = StellarRadius_SI;
				map_reference_properties["StellarRadius_norm"] = StellarRadius_norm;
				map_reference_properties["MassLossRate_SI"] = MassLossRate_SI;
				map_reference_properties["MassLossRate_norm"] = MassLossRate_norm;


			}


		}

		h5in.close();

	}

	// Now, we synchronise the map in the MPI-parallel case:
#ifdef parallel
	broadcast_data();
#endif

	// Compute magnetic field -- remark: currently, there are no derivatives in the phi-direction
	Br_ref.resize(Index::set(0,0,0),
			Index::set(mx_ref[0]+1, mx_ref[1], mx_ref[2]));
	Btheta_ref.resize(Index::set(0,0,0),
			Index::set(mx_ref[0], mx_ref[1]+1, mx_ref[2]));
	Bphi_ref.resize(Index::set(0,0,0),
			Index::set(mx_ref[0], mx_ref[1], mx_ref[2]+1));

	if(mpi_rank==0) {
		cout << " Computing magnetic field from vector-potential data" << endl;
	}

	// Compute Br from numerical curl
	for(int iphi=0; iphi<=mx_ref[2]; ++iphi) {
		for(int itheta=0; itheta<=mx_ref[1]; ++itheta) {
			double theta_R = thetaL_ref(itheta+1);
			double theta_L = thetaL_ref(itheta);
			for(int ir=0; ir<=mx_ref[0]+1; ++ir) {
				Br_ref(ir,itheta,iphi) = ((sin(theta_R)*Aphi_ref(ir,itheta+1,iphi) -
						sin(theta_L)*Aphi_ref(ir,itheta,iphi))/(theta_R - theta_L))/(radL(ir)*sin(thetaCen_ref(itheta)));
			}
		}
	}
	// Compute Btheta from numerical curl -- currently ignoring phi-dependecies
	for(int iphi=0; iphi<=mx_ref[2]; ++iphi) {
		for(int itheta=0; itheta<=mx_ref[1]+1; ++itheta) {
			for(int ir=0; ir<=mx_ref[0]; ++ir) {
				double rad_R = radL(ir+1);
				double rad_L = radL(ir);
				Btheta_ref(ir,itheta,iphi) = -((rad_R*Aphi_ref(ir+1,itheta,iphi) -
						rad_L*Aphi_ref(ir,itheta,iphi))/(rad_R - rad_L))/(radCen(ir));
			}
		}
	}

	Bphi_ref.clear();

	if(mpi_rank==0) {
		cout << " Remapping magnetic-field data " << endl;
	}

	NumMatrix<float,3> Br_old(Br_ref), Btheta_old(Btheta_ref), Bphi_old(Bphi_ref);
	// Remap magnetic-field to cell centers:
	for(int iphi=0; iphi<=mx_ref[2]; ++iphi) {
		for(int itheta=0; itheta<=mx_ref[1]; ++itheta) {
			for(int ir=0; ir<=mx_ref[0]; ++ir) {

				Br_ref(ir,itheta,iphi) = 0.5*(Br_old(ir,itheta,iphi) + Br_old(ir+1,itheta,iphi));
				Btheta_ref(ir,itheta,iphi) = 0.5*(Btheta_old(ir,itheta,iphi) + Btheta_old(ir,itheta+1,iphi));
				Bphi_ref(ir,itheta,iphi) = 0.5*(Bphi_old(ir,itheta,iphi) + Bphi_old(ir,itheta,iphi+1));

			}
		}
	}


	if(mpi_rank==0) {
		cout << " Remapping vector-potential data " << endl;
	}
	// Remap vector potential to cell centers:
	// This works using just one array, because we only overwrite the oldest values
	for(int iphi=0; iphi<=mx_ref[2]; ++iphi) {
		for(int itheta=0; itheta<=mx_ref[1]; ++itheta) {
			for(int ir=0; ir<=mx_ref[0]; ++ir) {

				Ar_ref(ir,itheta,iphi) = 0.25*(Ar_ref(ir,itheta,iphi) +
						Ar_ref(ir,itheta,iphi+1) + Ar_ref(ir,itheta+1,iphi) +
						Ar_ref(ir,itheta+1,iphi+1));
				dArdt_ref(ir,itheta,iphi) = 0.25*(dArdt_ref(ir,itheta,iphi) +
						dArdt_ref(ir,itheta,iphi+1) + dArdt_ref(ir,itheta+1,iphi) +
						dArdt_ref(ir,itheta+1,iphi+1));

				Atheta_ref(ir,itheta,iphi) = 0.25*(Atheta_ref(ir,itheta,iphi) +
						Atheta_ref(ir,itheta,iphi+1) + Atheta_ref(ir+1,itheta,iphi) +
						Atheta_ref(ir+1,itheta,iphi+1));
				dAthetadt_ref(ir,itheta,iphi) = 0.25*(dAthetadt_ref(ir,itheta,iphi) +
						dAthetadt_ref(ir,itheta,iphi+1) + dAthetadt_ref(ir+1,itheta,iphi) +
						dAthetadt_ref(ir+1,itheta,iphi+1));

				Aphi_ref(ir,itheta,iphi) = 0.25*(Aphi_ref(ir,itheta,iphi) +
						Aphi_ref(ir+1,itheta,iphi) + Aphi_ref(ir,itheta+1,iphi) +
						Aphi_ref(ir+1,itheta+1,iphi));
				dAphidt_ref(ir,itheta,iphi) = 0.25*(dAphidt_ref(ir,itheta,iphi) +
						dAphidt_ref(ir+1,itheta,iphi) + dAphidt_ref(ir,itheta+1,iphi) +
						dAphidt_ref(ir+1,itheta+1,iphi));

			}
		}
	}

	if(mpi_rank==0) {
		cout << " Done remapping " << endl;
	}

	//		cout << " After remap " << endl;
	//		for(int ir=0; ir<350; ir++) {
	//							cout << ir << " " << vr_ref(ir,10,0) << " " << vr_ref.getHigh(0) << endl;
	//						}



	return 0;
}


bool ReferenceSolution::check_property(std::string name_property) {
	if( map_reference_properties.count( name_property ) > 0 ) {
		return true;
	} else {
		return false;
	}
}

float ReferenceSolution::get_property(std::string name_property) {
//	std::map<char,int>::iterator it;
//	it = map_reference_properties.find(name_property);
	return map_reference_properties.find(name_property)->second;

}

int ReferenceSolution::get_gridIndexRad(double r_sph, bool use_rCen) {


	if(grid_type_rad==6) {
		double irVal = log(1.- (r_sph-grid_beg_rad)*(1.-grid_inc_rad)/grid_del_rad)/log(grid_inc_rad);
		if(use_rCen) {
			irVal -= 0.5;
		}
//		cout << " finding index " << r_sph << " " << grid_type_rad << " " << irVal <<  endl;
		i_refNearestL(0) = static_cast<int>(irVal);
	} else {
		int ir = i_refNearestL(0);
		if(ir<0) {
			ir=0;
		} else if (ir>mx_ref[0]) {
			ir = mx_ref[0];
		}
		double delta;

		if(use_rCen) {
			delta = r_sph - radCen[ir];
		} else {
			delta = r_sph - radL[ir];
		}

		if(delta > 0.) {
			while(delta > 0 && ir<mx_ref[0]) {
				ir++;
				if(use_rCen) {
					delta = r_sph - radCen[ir];
				} else {
					delta = r_sph - radL[ir];
				}
			}
			i_refNearestL(0) = ir-1;
		} else {
			while(delta < 0 && ir>0) {
				ir--;
				if(use_rCen) {
					delta = r_sph - radCen[ir];
				} else {
					delta = r_sph - radL[ir];
				}
			}
			i_refNearestL(0) = ir;
		}
	}



	// Capture cases, where point is outside of reference grid:
	if(i_refNearestL(0)<0) {
		i_refNearestL(0) = 0;
		// Set shifts so that ir=0 is used exclusively
		dist_fromRefL(0) = 0.;
		dist_refPoints(0) = 1.;
		dist_fromRefR(0) = dist_refPoints(0) - dist_fromRefL(0);
	} else if (i_refNearestL(0) > mx_ref[0]-1) {
		i_refNearestL(0) = mx_ref(0)-1;
		// set shifts so that ir=mx_ref[0] is used exclusively
		dist_refPoints(0) = 1.;
		dist_fromRefR(0) = 0.;
		dist_fromRefL(0) = dist_refPoints(0) - dist_fromRefR(0);
	} else {
		// general case ->  i_refNearestL can be used as it is
		float rad_refNearestL = (use_rCen) ? radCen[i_refNearestL(0)] : radL[i_refNearestL(0)]; // radial position of nearest reference point
		dist_fromRefL(0) = r_sph - rad_refNearestL;
		dist_refPoints(0) = (use_rCen) ? radCen[i_refNearestL(0)+1] : radL[i_refNearestL(0)+1];
		dist_refPoints(0) -= rad_refNearestL;
		dist_fromRefR(0) = dist_refPoints(0) - dist_fromRefL(0);
	}

//	cout << " final index " << i_refNearestL(0) << endl;

	return i_refNearestL(0);
}


int ReferenceSolution::get_gridIndexTheta(double thetaPos, bool use_thetaCen) {

	if(grid_type_theta==0) {
		double itVal = (thetaPos-grid_beg_theta)/grid_del_theta;
		if(use_thetaCen) {
			itVal -= 0.5;
		}
		i_refNearestL(1) = static_cast<int>(itVal);

	} else {

		// find position iteratively
		int it = i_refNearestL(1);
		double delta;

		if(use_thetaCen) {
			delta = thetaPos - thetaCen_ref[it];
		} else {
			delta = thetaPos - thetaL_ref[it];
		}

		if(delta > 0.) {
			while(delta > 0 && it<mx_ref[1]) {
				it++;
				if(use_thetaCen) {
					delta = thetaPos - thetaCen_ref[it];
				} else {
					delta = thetaPos - thetaL_ref[it];
				}
			}
			i_refNearestL(1) = it-1;
		} else {
			while(delta < 0 && it>0) {
				it--;
				if(use_thetaCen) {
					delta = thetaPos - thetaCen_ref[it];
				} else {
					delta = thetaPos - thetaL_ref[it];
				}

			}
			i_refNearestL(1) = it;
		}
	}



	// Capture cases, where point is outside of reference grid:
	if(i_refNearestL(1)<0) {
		i_refNearestL(1) = 0;
		// Set shifts so that ir=0 is used exclusively
		dist_fromRefL(1) = 0.;
		dist_refPoints(1) = 1.;
		dist_fromRefR(1) = dist_refPoints(1) - dist_fromRefL(1);

	} else if (i_refNearestL(1) == mx_ref[1]) {

		// If point at boundary - cannot use that one, because non exists above that
		i_refNearestL(1) = mx_ref[1] - 1;
		float refPos = (use_thetaCen) ? thetaCen_ref[i_refNearestL(1)] : thetaL_ref[i_refNearestL(1)];
		dist_fromRefL(1) = thetaPos - refPos;
		dist_refPoints(1) = (use_thetaCen) ? thetaCen_ref[i_refNearestL(1)+1] : thetaL_ref[i_refNearestL(1)+1];
		dist_refPoints(1) -= refPos;
		dist_fromRefR(1) = dist_refPoints(1) - dist_fromRefL(1);

	} else if (i_refNearestL(1)>mx_ref[1]) {

		i_refNearestL(1) = mx_ref[1]-1;
		dist_refPoints(1) = 1.;
		dist_fromRefR(1) = 0.;
		dist_fromRefL(0) = dist_refPoints(1) - dist_fromRefR(1);

	} else {

		// general case ->  i_refNearestL can be used as it is
		float theta_refNearestL = (use_thetaCen) ? thetaCen_ref[i_refNearestL(1)] : thetaL_ref[i_refNearestL(1)];
		dist_fromRefL(1) = thetaPos - theta_refNearestL;
		dist_refPoints(1) = (use_thetaCen) ? thetaCen_ref[i_refNearestL(1)+1] : thetaL_ref[i_refNearestL(1)+1];
		dist_refPoints(1) -= theta_refNearestL;
		//		refDels(1) = thetaRef[refPoints(1)+1] - thetaRef[refPoints(1)];
		dist_fromRefR(1) = dist_refPoints(1) - dist_fromRefL(1);
	}

	return i_refNearestL(1);
}

double ReferenceSolution::get_interpolated(NumMatrix<float, 3> &field, int i_rad, int j_theta, int k_phi) {
	// Ignore phi dependence until further notice
	double phi = (field(i_rad,j_theta,k_phi)*dist_fromRefR(0)*dist_fromRefR(1) +
			field(i_rad  ,j_theta+1,k_phi)*dist_fromRefR(0)*dist_fromRefL(1) +
			field(i_rad+1,j_theta  ,k_phi)*dist_fromRefL(0)*dist_fromRefR(1) +
			field(i_rad+1,j_theta+1,k_phi)*dist_fromRefL(0)*dist_fromRefL(1))/(dist_refPoints(0)*dist_refPoints(1));

//	cout << " values ";
//	cout << field(i_rad,j_theta,k_phi) << " ";
//	cout << field(i_rad+1,j_theta,k_phi) << " ";
//	cout << field(i_rad,j_theta+1,k_phi) << " ";
//	cout << field(i_rad+1,j_theta+1,k_phi) << " ";
//	cout << endl;
//	cout << "        ";
//	cout << dist_fromRefR(0) << " " << dist_fromRefR(1) << " " << dist_fromRefL(0) << " " << dist_fromRefL(1) << endl;


	return phi;
}


NumArray<double> ReferenceSolution::transformVec_SphToCart(double v_r, double v_theta, double v_phi, double theta, double phi) {
	NumArray<double> cartVector(3);

	double sinTheta = sin(theta);
	double cosTheta = cos(theta);
	double sinPhi = sin(phi);
	double cosPhi = cos(phi);

	cartVector(0) =  v_r*sinTheta*cosPhi + v_theta*cosTheta*cosPhi - v_phi*sinPhi;
	cartVector(1) =  v_r*sinTheta*sinPhi + v_theta*cosTheta*sinPhi + v_phi*cosPhi;
	cartVector(2) = v_r*cosTheta - v_theta*sinTheta;

	return cartVector;
}


void ReferenceSolution::tilt_solution(NumArray<double> &vector, double direction) {
	tilt_solution(vector(0), vector(1), vector(2), direction);
}

void ReferenceSolution::tilt_solution(double &vec_x, double &vec_y, double &vec_z, double direction) {
	double theta_tilt_curr = theta_tilt*direction;

	double vec_x_old = vec_x;
	double vec_z_old = vec_z;

	vec_x =  vec_x_old*cos(theta_tilt_curr) + vec_z_old*sin(theta_tilt_curr);
	vec_z = -vec_x_old*sin(theta_tilt_curr) + vec_z_old*cos(theta_tilt_curr);
}


void ReferenceSolution::get_ref_hydro(NumArray<double> &refSolution, double x_rel, double y_rel, double z_rel) {
	assert(refSolution.getLength()>4);

	// tilt Cartesian position vector if necessary:
	if(this->do_tilt) {
		tilt_solution(x_rel, y_rel, z_rel);
	}


	double r_sph = sqrt(sqr(x_rel) + sqr(y_rel) + sqr(z_rel));
	double theta = atan2 (sqrt(sqr(x_rel) + sqr(y_rel)), z_rel);
	double phi = atan2(y_rel, x_rel);

	int i_rad = get_gridIndexRad(r_sph, 1);
	int j_theta = get_gridIndexTheta(theta, 1);

	// Get density
	refSolution(0) = get_interpolated(rho_ref, i_rad, j_theta, 0);

//	cout << " Bef " << endl;
//			for(int ir=0; ir<350; ir++) {
//					cout << ir << " " << vr_ref(ir,10,0) << " " << vr_ref.getHigh(0) << endl;
//				}


	// Get spherical velocity components
	double v_rad = get_interpolated(vr_ref, i_rad, j_theta, 0);
	double v_theta = get_interpolated(vtheta_ref, i_rad, j_theta, 0);
	double v_phi = get_interpolated(vphi_ref, i_rad, j_theta, 0);

//	cout << " velo ";
//	cout << v_rad << " " << v_theta << " " << v_phi << endl;
//	cout << " ";
//	cout << i_rad << " " << j_theta << " " << vr_ref(i_rad,j_theta,0) << " " << endl;


//	for(int ir=0; ir<350; ir++) {
//		cout << ir << " " << vr_ref(ir,10,0) << " " << vr_ref.getHigh(0) << endl;
//	}
//	exit(3);

	NumArray<double> cartVector(3);
	cartVector = transformVec_SphToCart(v_rad, v_theta, v_phi, theta, phi);

	refSolution(1) = cartVector(0);
	refSolution(2) = cartVector(1);
	refSolution(3) = cartVector(2);

//	cout << " ref " << x_rel << " " << y_rel << " " << z_rel << " ";
//	cout << refSolution(1) << " ";
//	cout << refSolution(2) << " ";
//	cout << refSolution(3) << " ";
//	cout << endl;

	// need to tilt velocity vector, too
	if(this->do_tilt) {
		tilt_solution(refSolution(1), refSolution(2), refSolution(3),-1.);
	}

	// Setting constant temperature (unfortunately, we do not know about the normalisation...)
	refSolution(4) = 1.;





}

void ReferenceSolution::get_ref_magnetic(NumArray<double> &refMag, NumArray<double> &refPot,
		NumArray<double> &d_refPot_dt, double x_rel, double y_rel, double z_rel) {
	assert(refMag.getLength()>2);
	assert(refPot.getLength()>2);
	assert(d_refPot_dt.getLength()>2);

	// tilt Cartesian position vector if necessary:
	if(this->do_tilt) {
		tilt_solution(x_rel, y_rel, z_rel);
	}


	double r_sph = sqrt(sqr(x_rel) + sqr(y_rel) + sqr(z_rel));
	double theta = atan2 (sqrt(sqr(x_rel) + sqr(y_rel)), z_rel);
	double phi = atan2(y_rel, x_rel);

	int i_rad = get_gridIndexRad(r_sph, 1);
	int j_theta = get_gridIndexTheta(theta, 1);


	// Get magnetic field components
	double B_rad = get_interpolated(Br_ref, i_rad, j_theta, 0);
	double B_theta = get_interpolated(Btheta_ref, i_rad, j_theta, 0);
	double B_phi = get_interpolated(Bphi_ref, i_rad, j_theta, 0);

	NumArray<double> cartVector(3);
	cartVector = transformVec_SphToCart(B_rad, B_theta, B_phi, theta, phi);

	refMag(0) = cartVector(0);
	refMag(1) = cartVector(1);
	refMag(2) = cartVector(2);

	// need to tilt velocity vector, too
	if(this->do_tilt) {
		tilt_solution(refMag(0), refMag(1), refMag(2),-1.);
	}



	// Get vector-potential components
	double A_rad = get_interpolated(Ar_ref, i_rad, j_theta, 0);
	double A_theta = get_interpolated(Atheta_ref, i_rad, j_theta, 0);
	double A_phi = get_interpolated(Aphi_ref, i_rad, j_theta, 0);

	cartVector = transformVec_SphToCart(A_rad, A_theta, A_phi, theta, phi);

	refPot(0) = cartVector(0);
	refPot(1) = cartVector(1);
	refPot(2) = cartVector(2);

	if(this->do_tilt) {
		tilt_solution(refPot(0), refPot(1), refPot(2),-1.);
	}


	// get time-rate of change of vector-potential components
	double dA_rad_dt = get_interpolated(dArdt_ref, i_rad, j_theta, 0);
	double dA_theta_dt = get_interpolated(dAthetadt_ref, i_rad, j_theta, 0);
	double dA_phi_dt = get_interpolated(dAphidt_ref, i_rad, j_theta, 0);

	cartVector = transformVec_SphToCart(dA_rad_dt, dA_theta_dt, dA_phi_dt, theta, phi);


	d_refPot_dt(0) = cartVector(0);
	d_refPot_dt(1) = cartVector(1);
	d_refPot_dt(2) = cartVector(2);

	if(this->do_tilt) {
		tilt_solution(d_refPot_dt(1), d_refPot_dt(2), d_refPot_dt(3),-1.);
	}

}



