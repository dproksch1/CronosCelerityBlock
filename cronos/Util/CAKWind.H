#ifndef CAK_WIND_H
#define CAK_WIND_H

#include "data.H"
#include "normalisation.H"

class CAKWind {
public:
	CAKWind(normalisation &norm,
			const Quantity &LStar_SI, const Quantity &tempStar, const Quantity &rad_star_SI, double k_star,
			double alpha_star, bool strong_coupling);
//	CAKWind(normalisation &norm, std::vector<Quantity> posStar_SI,
//			const Quantity &LStar_SI, const Quantity &tempStar, const Quantity &rad_star_SI, double k_star,
//			double alpha_star, bool strong_coupling);
	/**Compute CAK acceleration just using CAK parameters for the star
	 * */
	double get_AccelAbs(Data &gdata, NumArray<double> &posStar, int i, int j, int k);

	/**Compute CAK acceleration with arbitrary CAK parameters, e.g., from the stellar wind */
	double get_AccelAbs(Data &gdata, NumArray<double> &posStar, int i, int j, int k, double k_wind, double alpha_wind);

	/**Compute CAK acceleration just using CAK parameters for the star
	 * */
	double get_AccelVec(Data &gdata, NumArray<double> &CAK_Accel, NumArray<double> &posStar, int i, int j, int k);

	/**Compute CAK acceleration with arbitrary CAK parameters, e.g., from the stellar wind */
	double get_AccelVec(Data &gdata, NumArray<double> &CAK_Accel, NumArray<double> &posStar, int i, int j, int k, double k_wind, double alpha_wind);

	/**Compute continuum acceleration*/
	double get_contAccel(Data &gdata, NumArray<double> &posStar, int ix, int iy, int iz);
private:
	/**Compute radial gradient of abs value of velocity*/
	double get_radialVelocityGradient(Data &gdata, int i, int j, int k);

	/**Compute the finite-disk correction factor - vectorial result */
	double get_FDFactor(Data &gdata, NumArray<double> &fd_factor, double alphaChoice, int ix, int iy, int iz);

	/**Get weight for Simpson integration*/
	double get_weightSimpson(int i_Simpson, int n_Simpson);

	/** Position of star in normalised units */
	NumArray<double> posStar_num;
	/**Properties of star (luminosity, effective temperature, and radius)*/
	double LStar_num, TeffStar_num, radStar_num;
	double kStar, alphaStar;
	double gRad_pre_num;//, gRad_wind_num;
	double gRad_continuum_num;
	double t_exp_num;
	double w_Simpson_bound, w_Simpson_cen, w_Simpson_inter;
	bool has_strongCoupling;
};


#endif
