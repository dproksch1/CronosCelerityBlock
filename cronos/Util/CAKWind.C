#include "CAKWind.H"

CAKWind::CAKWind(normalisation &norm,
		const Quantity &LStar_SI, const Quantity &TeffStar_SI, const Quantity &rad_star_SI,
		double k_star, double alpha_star, bool strong_coupling) {

//	for (unsigned i_dir=0; i_dir < posStar_SI.size(); i_dir++) {
//		posStar_num[i_dir] = norm.get_num(posStar_SI[i_dir]);
//	}

	this->LStar_num = norm.get_num(LStar_SI);
	this->TeffStar_num = norm.get_num(TeffStar_SI);
	this->radStar_num = norm.get_num(rad_star_SI);

	this->kStar = k_star;

	this->alphaStar = alpha_star;
	this->has_strongCoupling = strong_coupling;

	// Helium fraction in wind -- currently assumed to be constant
	double frac_He = 0.1;
	double I_He; // Ionisation of helium
	if(TeffStar_SI > 30000) {
		I_He = 2.;
	} else {
		I_He = 1.;
	}

	Quantity sigma_Th = CRONOS_CONSTANTS::ThomsonCrossSection;
	Quantity m_H = CRONOS_CONSTANTS::HydrogenMass;
	Quantity kB = CRONOS_CONSTANTS::BoltzmannConstant;
	Quantity c_light = CRONOS_CONSTANTS::SpeedOfLight;

	// Compute specific electron opacity
	Quantity sigma_e = sigma_Th/m_H*(1.+I_He*frac_He)/(1.+4.*frac_He);
	// Normalise as specfic cross section
	double sigma_e_num = norm.get_num(norm.CROSS_SEC_S, sigma_e);

	// Now effects of normalisation need to be considered

	// thermal speed
//	Quantity v_therm = sqrt(2*kB*TeffStar/m_H)*CRONOS_BASE_UNITS::Meter/CRONOS_BASE_UNITS::Second;
	Quantity v_therm = sqrt(2*kB*TeffStar_SI/m_H);
	double v_therm_num = norm.get_num(v_therm);


	Quantity gRad_pre_SI = LStar_SI/(4.*pi*sqr(rad_star_SI))*sigma_e/c_light;
//	Quantity gRad_wind_phys = kWind*luminosity_star/(4.*pi*sqr(rad_star_SI))*sigma_e/c_light *
//			CRONOS_BASE_UNITS::Meter/sqr(CRONOS_BASE_UNITS::Second);

	// Store numerical pre-factor for t (see CAK 1975 Eq. (5)); Take care of additional factor 2 by Klaus(?)
	t_exp_num = sigma_e_num*v_therm_num;


	// Get numerical pre-factor for acceleration:
	gRad_pre_num= norm.get_num( gRad_pre_SI );
//	gRad_star_num = norm.get_num(norm.ACCEL, gRad_star_phys);
//	gRad_wind_num = norm.get_num(norm.ACCEL, gRad_wind_phys);

	Quantity gRad_continuum_SI =  LStar_SI*sigma_e/(4*pi*c_light);

	// Get numerical pre-factor for continuum radiation
	gRad_continuum_num = norm.get_num( gRad_continuum_SI );

//	cout << " Continuum radiation factor " << gRad_continuum_SI << " " << gRad_continuum_num << endl;
//	exit(39);


	w_Simpson_bound = 1./6.;
	w_Simpson_cen = 4./6.;
	w_Simpson_inter = 2./6.;

//	gRad_star_phys *= pow(1./(sigma_e*rho_norm*v_therm)*velo_phys/len_phys, alphaStar);

//	cerr << " Careful: continuum radiation should be absorbed in gravity term " << endl;
//	exit(3);

}


double CAKWind::get_contAccel(Data &gdata, NumArray<double> &posStar, int ix, int iy, int iz) {
	// Compute distance from star
	double dist = sqrt(sqr(gdata.getCen_x(ix) - posStar_num[0])
			+ sqr(gdata.getCen_y(iy) - posStar_num[1])
			+ sqr(gdata.getCen_z(iz) - posStar_num[2]));
	double continuum_force = gRad_continuum_num/sqr(dist);
	return continuum_force;
}

double CAKWind::get_AccelAbs(Data &gdata, NumArray<double> &posStar, int ix, int iy, int iz) {
	return get_AccelAbs(gdata, posStar, ix, iy, iz, kStar, alphaStar);
}


double CAKWind::get_AccelAbs(Data &gdata, NumArray<double> &posStar, int ix, int iy, int iz, double kChoice, double alphaChoice) {

//	posStar_num = posStar;
//
//	// Compute distance from star
//	double dist = sqrt(sqr(gdata.getCen_x(ix) - posStar_num[0])
//			+ sqr(gdata.getCen_y(iy) - posStar_num[1])
//			+ sqr(gdata.getCen_z(iz) - posStar_num[2]));
//
////	 gRad_w1_s1_norm -> K_CAK_Star
////			 gRad_w2_s1_norm -> K_CAK_Wind
//	if(!has_strongCoupling) {
//		kChoice = kStar;
//		alphaChoice = alphaStar;
//	}
//
//	double baseForce = gRad_pre_num*sqr(radStar_num/dist);
//
////	double optical_depth = t_exp_num*gdata.om[0](ix,iy,iz);
//	double optical_depth = t_exp_num*gdata.om[0](ix,iy,iz)/std::abs(get_radialVelocityGradient(gdata,ix,iy,iz));
//
//	double forceMultiplier = kChoice * pow(optical_depth, -alphaChoice);
//
//	double CAKForce = baseForce*forceMultiplier;
//
//
//	bool use_FDCorrection = true;
//
//	NumArray<double> FDCorrection_Factor_vec(3);
//	FDCorrection_Factor_vec.clear();
//
//	// Apply finite-disk correction factor
//	if(use_FDCorrection) {
//		double FDCorrection_Factor =  get_FDFactor(gdata, FDCorrection_Factor_vec, alphaChoice, ix, iy, iz);
////		cout << " comp " << FDCorrection_Factor << " ";
////		cout << sqrt(sqr(FDCorrection_Factor_vec(0)) + sqr(FDCorrection_Factor_vec(1)) + sqr(FDCorrection_Factor_vec(2)));
////		cout << endl;
//		CAKForce *= FDCorrection_Factor;
//	}
	NumArray<double> CAK_Accel(3);
	double CAKForce = get_AccelVec(gdata, CAK_Accel, posStar, ix, iy, iz, kChoice, alphaChoice);

	return CAKForce;
}


double CAKWind::get_AccelVec(Data &gdata, NumArray<double> &CAK_Accel, NumArray<double> &posStar, int ix, int iy, int iz) {
	return get_AccelVec(gdata, CAK_Accel, posStar, ix, iy, iz, kStar, alphaStar);
}


double CAKWind::get_AccelVec(Data &gdata, NumArray<double> &CAKForce_vec, NumArray<double> &posStar, int ix, int iy, int iz, double kChoice, double alphaChoice) {

	posStar_num = posStar;

	// Compute distance from star
	double dist = sqrt(sqr(gdata.getCen_x(ix) - posStar_num[0])
			+ sqr(gdata.getCen_y(iy) - posStar_num[1])
			+ sqr(gdata.getCen_z(iz) - posStar_num[2]));

//	 gRad_w1_s1_norm -> K_CAK_Star
//			 gRad_w2_s1_norm -> K_CAK_Wind
	if(!has_strongCoupling) {
		kChoice = kStar;
		alphaChoice = alphaStar;
	}

	double baseForce = gRad_pre_num*sqr(radStar_num/dist);

//	double optical_depth = t_exp_num*gdata.om[0](ix,iy,iz);
	double optical_depth = t_exp_num*gdata.om[0](ix,iy,iz)/std::abs(get_radialVelocityGradient(gdata,ix,iy,iz));

	double forceMultiplier = kChoice * pow(optical_depth, -alphaChoice);

	double CAKForce = baseForce*forceMultiplier;


	bool use_FDCorrection = true;

	NumArray<double> FDCorrection_Factor_vec(3);
	FDCorrection_Factor_vec.clear();

	// Apply finite-disk correction factor
	if(use_FDCorrection) {
		double FDCorrection_Factor =  get_FDFactor(gdata, FDCorrection_Factor_vec, alphaChoice, ix, iy, iz);
//		cout << " comp " << FDCorrection_Factor << " ";
//		cout << sqrt(sqr(FDCorrection_Factor_vec(0)) + sqr(FDCorrection_Factor_vec(1)) + sqr(FDCorrection_Factor_vec(2)));
//		cout << endl;
		CAKForce_vec = FDCorrection_Factor_vec*CAKForce;
		CAKForce *= FDCorrection_Factor;
//		cout << " vec: ";
//		cout << CAKForce_vec(0) << " ";
//		cout << CAKForce_vec(1) << " ";
//		cout << CAKForce_vec(2) << " ";
//		cout << endl;
	} else {
		// Without finite disk correction factor, the force points away from the star
		CAKForce_vec(0) = gdata.getCen_x(ix) - posStar_num(0);
		CAKForce_vec(1) = gdata.getCen_y(iy) - posStar_num(1);
		CAKForce_vec(2) = gdata.getCen_z(iz) - posStar_num(2);
		CAKForce_vec *= CAKForce/dist;
	}

	return CAKForce;
}


double CAKWind::get_FDFactor(Data &gdata, NumArray<double> &fd_factor,
		double alphaChoice, int ix, int iy, int iz) {

	// Compute distance from star
	double dist = sqrt(sqr(gdata.getCen_x(ix) - posStar_num[0])
			+ sqr(gdata.getCen_y(iy) - posStar_num[1])
			+ sqr(gdata.getCen_z(iz) - posStar_num[2]));
	// Compute size of stellar disk
	double mumin = sqrt(1-sqr(radStar_num/dist));

	// Factor in front of finite-disk correction factor (1/(1-mu_Star**2))*1/pi = r**2/RStar**2 * 1/pi
	double constInteg = sqr(dist/radStar_num)/M_PI;


	int q_sx(1), q_sy(2), q_sz(3);

	double inv_del_x = 1./(gdata.getCen_x(ix+1)-gdata.getCen_x(ix-1));
	double inv_del_y = 1./(gdata.getCen_y(iy+1)-gdata.getCen_y(iy-1));
	double inv_del_z = 1./(gdata.getCen_z(iz+1)-gdata.getCen_z(iz-1));

	REAL dv_xdx= (gdata.om[q_sx](ix+1,iy,iz)-gdata.om[q_sx](ix-1,iy,iz)) * inv_del_x;
	REAL dv_ydx= (gdata.om[q_sy](ix+1,iy,iz)-gdata.om[q_sy](ix-1,iy,iz)) * inv_del_x;
	REAL dv_zdx= (gdata.om[q_sz](ix+1,iy,iz)-gdata.om[q_sz](ix-1,iy,iz)) * inv_del_x;
	REAL dv_xdy= (gdata.om[q_sx](ix,iy+1,iz)-gdata.om[q_sx](ix,iy-1,iz)) * inv_del_y;
	REAL dv_ydy= (gdata.om[q_sy](ix,iy+1,iz)-gdata.om[q_sy](ix,iy-1,iz)) * inv_del_y;
	REAL dv_zdy= (gdata.om[q_sz](ix,iy+1,iz)-gdata.om[q_sz](ix,iy-1,iz)) * inv_del_y;
	REAL dv_xdz= (gdata.om[q_sx](ix,iy,iz+1)-gdata.om[q_sx](ix,iy,iz-1)) * inv_del_z;
	REAL dv_ydz= (gdata.om[q_sy](ix,iy,iz+1)-gdata.om[q_sy](ix,iy,iz-1)) * inv_del_z;
	REAL dv_zdz= (gdata.om[q_sz](ix,iy,iz+1)-gdata.om[q_sz](ix,iy,iz-1)) * inv_del_z;

	// Get distances in respective plains
	REAL dist_xy = sqrt(sqr(gdata.getCen_x(ix) - posStar_num[0])
			+ sqr(gdata.getCen_y(iy) - posStar_num[1]));
	REAL dist_yz = sqrt(sqr(gdata.getCen_z(iz) - posStar_num[2])
			+ sqr(gdata.getCen_y(iy) - posStar_num[1]));
	REAL dist_xz = sqrt(sqr(gdata.getCen_x(ix) - posStar_num[0])
			+ sqr(gdata.getCen_z(iz) - posStar_num[2]));


	// Get local unit vectors
	// Compute normal vector from centre of solar surface
	NumArray<double> normStar(3), perpStar(3), binormStar(3);
	normStar[0] = (posStar_num[0]-gdata.getCen_x(ix))/dist;
	normStar[1] = (posStar_num[1]-gdata.getCen_y(iy))/dist;
	normStar[2] = (posStar_num[2]-gdata.getCen_z(iz))/dist;

	// define perpendicular vector
	if((dist_xy >= dist_xz) && (dist_xy >= dist_yz)){
		perpStar[0] = (posStar_num[1]-gdata.getCen_y(iy))/dist_xy;
		perpStar[1] = -(posStar_num[0]-gdata.getCen_x(ix))/dist_xy;
		perpStar[2] = 0;
	} else if (dist_xz >= dist_yz){
		perpStar[0] = -(posStar_num[2]-gdata.getCen_z(iz))/dist_xz;
		perpStar[1] =  0;
		perpStar[2] =  (posStar_num[0]-gdata.getCen_x(ix))/dist_xz;
	} else {
		perpStar[0] = 0;
		perpStar[1] =  (posStar_num[2]-gdata.getCen_z(iz))/dist_yz;
		perpStar[2] = -(posStar_num[1]-gdata.getCen_y(iy))/dist_yz;
	}

	// Compute vector perpendicular on both perp and N:
	binormStar[0] = normStar[1]*perpStar[2] - normStar[2]*perpStar[1];
	binormStar[1] = normStar[2]*perpStar[0] - normStar[0]*perpStar[2];
	binormStar[2] = normStar[0]*perpStar[1] - normStar[1]*perpStar[0];



	double radialGradient = get_radialVelocityGradient(gdata, ix, iy, iz);


	// Integrate over stellar surface:
	int mu_steps(2), phi_steps(4);
	int Simpson_steps = 2*mu_steps + 1;
	double g_rad = 0.;
	REAL integral = 0.;
	REAL del_mu = (1.-mumin)/(1.*mu_steps);
	REAL del_mu_Simpson = (1.-mumin)/(Simpson_steps-1.);
	REAL del_phi = 2.*M_PI/(1.*phi_steps);
	NumArray<double> perpLoc(3), normLoc(3);
	fd_factor.clear();

	// Hardcoded to four steps in phi-direction with different respective unit vectors
	for(int iPhi=0; iPhi<phi_steps; ++iPhi) {

		double phi = del_phi*iPhi; // + 0.5*iPhi;

		// Compute a vector perpendicular to normal vector for all directions at specific phi-direction
		for(int i_dir=0; i_dir<3; ++i_dir) {
			perpLoc(i_dir) = perpStar(i_dir)*cos(phi) + binormStar(i_dir)*sin(phi);
//			perpLoc(0) = perpStar(0)*cos(phi) + binormStar(0)*sin(phi);
		}


		REAL mu = mumin;

		for(int i_mu=0; i_mu<mu_steps; ++i_mu) {

		}

		// Integrate using Simpson's rule
		for(int i_simpsonStep=0; i_simpsonStep<Simpson_steps; ++i_simpsonStep) {

			double mu = mumin + del_mu_Simpson*i_simpsonStep;
			if (sqrt(sqr((1-sqr(mu)))) < 1e-14) {mu=1;}

			// Normal vector to local patch from normal vector to center of stellar disc
			// (normStar) and perpendicular vector on stellar surface (perpLoc)
			normLoc(0) = (sqrt(1-sqr(mu))*perpLoc(0)+mu*normStar[0]);
			normLoc(1) = (sqrt(1-sqr(mu))*perpLoc(1)+mu*normStar[1]);
			normLoc(2) = (sqrt(1-sqr(mu))*perpLoc(2)+mu*normStar[2]);

			double projected_grad = sqr(normLoc[0])*dv_xdx + normLoc[0]*normLoc[1]*(dv_xdy+dv_ydx) +
					sqr(normLoc[1])*dv_ydy + normLoc[0]*normLoc[2]*(dv_xdz+dv_zdx) +
					sqr(normLoc[2])*dv_zdz + normLoc[1]*normLoc[2]*(dv_ydz+dv_zdy);
			projected_grad = std::abs(projected_grad); // Use absolute value

			double integrand = pow(std::abs(projected_grad/radialGradient),alphaChoice);

			double w_Simpson = get_weightSimpson(i_simpsonStep, Simpson_steps);

			integral += integrand*w_Simpson*mu;

			fd_factor(0) += -integrand*w_Simpson*normLoc(0);
			fd_factor(1) += -integrand*w_Simpson*normLoc(1);
			fd_factor(2) += -integrand*w_Simpson*normLoc(2);

		}

	}
	integral *= del_phi*del_mu*constInteg;
	fd_factor *= del_phi*del_mu*constInteg;

	return integral;


}




double CAKWind::get_weightSimpson(int i_Simpson, int n_Simpson) {
	if(i_Simpson==0 || i_Simpson==n_Simpson-1) {
		return w_Simpson_bound;
	} else if (i_Simpson%2==1) {
		return w_Simpson_cen;
	} else {
		return w_Simpson_inter;
	}
}


double CAKWind::get_radialVelocityGradient(Data &gdata, int i, int j, int k) {

	NumArray<double> dist(3);

	int q_sx(1), q_sy(2), q_sz(3);

	dist(0) = gdata.getCen_x(i) - posStar_num[0];
	dist(1) = gdata.getCen_y(j) - posStar_num[1];
	dist(2) = gdata.getCen_z(k) - posStar_num[2];

	double distAbs = sqrt(sqr(dist(0)) + sqr(dist(1)) + sqr(dist(2)));

	NumArray<double> grad_u(3);

	grad_u(0) = (sqrt(sqr(gdata.om[q_sx](i+1,j,k)) + sqr(gdata.om[q_sy](i+1,j,k)) + sqr(gdata.om[q_sz](i+1,j,k))) -
			sqrt(sqr(gdata.om[q_sx](i-1,j,k)) + sqr(gdata.om[q_sy](i-1,j,k)) + sqr(gdata.om[q_sz](i-1,j,k))))
			/ (gdata.getCen_x(i+1)-gdata.getCen_x(i-1));
	grad_u(1) = (sqrt(sqr(gdata.om[q_sx](i,j+1,k)) + sqr(gdata.om[q_sy](i,j+1,k)) + sqr(gdata.om[q_sz](i,j+1,k))) -
			sqrt(sqr(gdata.om[q_sx](i,j-1,k)) + sqr(gdata.om[q_sy](i,j-1,k)) + sqr(gdata.om[q_sz](i,j-1,k))))
			/ (gdata.getCen_y(j+1)-gdata.getCen_y(j-1));
	grad_u(2) = (sqrt(sqr(gdata.om[q_sx](i,j,k+1)) + sqr(gdata.om[q_sy](i,j,k+1)) + sqr(gdata.om[q_sz](i,j,k+1))) -
			sqrt(sqr(gdata.om[q_sx](i,j,k-1)) + sqr(gdata.om[q_sy](i,j,k-1)) + sqr(gdata.om[q_sz](i,j,k-1))))
			/ (gdata.getCen_z(k+1)-gdata.getCen_z(k-1));

	double grad_rad = (grad_u(0)*dist(0) +  grad_u(1)*dist(1) + grad_u(2)*dist(2))/distAbs;

	return grad_rad;
}

