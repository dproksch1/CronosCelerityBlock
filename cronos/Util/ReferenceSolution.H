#ifndef CRONOS_REFERENCE_SOLUTION_H
#define CRONOS_REFERENCE_SOLUTION_H

#include "matrix.H"
#include <string>
#include <map>

class ReferenceSolution {
public:
	/**Constructor for reference-solution class - need to give full path to solution file*/
	ReferenceSolution(std::string full_path, int rank=0);
	/**Check if property is available*/
	bool check_property(std::string name_property);
	/**Get value for specific property*/
	float get_property(std::string name_property);
	/**Get reference value of standard variables at cell centre*/
	void get_ref_hydro(NumArray<double> &refSolution, double x_rel, double y_rel, double z_rel);
	/**Get reference value of of vector potential */
	void get_ref_magnetic(NumArray<double> &refMag, NumArray<double> &refPot, NumArray<double> &d_refPot_dt, double x_rel, double y_rel, double z_rel);
	/**Selt tilt angle for solution in theta direction*/
	void set_tiltAngle(double tilt_angle);
private:
	/**Test if file exists*/
	bool file_exists(const std::string& filename);

	/**Load reference data from hdf5 file*/
	int load_file();

	/**Send loaded data from rank 0 to other MPI ranks*/
	void broadcast_data();

	/**Transfer map entry from rank 0 to other MPI ranks*/
	void broadcast_map_entry(std::string &name_property, float &value_property);

	/**get nearest left-handed index on radial grid*/
	int get_gridIndexRad(double r_sph, bool use_rCen);

	/**get nearest right-handed index on theta grid*/
	int get_gridIndexTheta(double theta, bool use_thetaCen);

	/**Compute interpolated point value*/
	double get_interpolated(NumMatrix<float, 3> &field, int i_rad, int j_theta, int k_phi);

	/**tilt solution by given angle in theta*/
	void tilt_solution(NumArray<double> &vector, double direction);

	void tilt_solution(double &vec_x, double &vec_y, double &vec_z, double direction=1.);

	/**Transform vector from spherical to Cartesian coordinates*/
	NumArray<double> transformVec_SphToCart(double v_r, double v_theta, double v_phi, double thata, double phi);

	NumMatrix<float, 3> rho_ref, vr_ref, vtheta_ref, vphi_ref;
	NumMatrix<float, 3> Ar_ref, Atheta_ref, Aphi_ref;
	NumMatrix<float, 3> Br_ref, Btheta_ref, Bphi_ref;
	NumMatrix<float, 3> dArdt_ref, dAthetadt_ref, dAphidt_ref;

	NumArray<float> radCen, radL, thetaCen_ref, thetaL_ref;
	/**Variables for local interpolation*/
	NumArray<float> dist_refPoints, dist_fromRefL, dist_fromRefR;

	/**Storage for nearest grid point left of point of interest*/
	NumArray<int> i_refNearestL;

	NumArray<int> mx_ref;

	std::map<std::string,float> map_reference_properties;

	std::string full_path;

	/**tilt angle if tilted solution is desired*/
	double theta_tilt;

	float grid_beg_rad, grid_del_rad, grid_inc_rad;
	float grid_beg_theta, grid_end_theta, grid_del_theta;

	int grid_type_rad, grid_type_theta;

	int i_rad_store;

	int mpi_rank;

	bool do_tilt;
};


#endif
