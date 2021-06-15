#ifndef CRONOS_ALFVEN_WINGS_SPHERICAL_H
#define CRONOS_ALFVEN_WINGS_SPHERICAL_H

#include "problem.H"

// using namespace std;

/*
  Alle Benutzerklassen erben einfach von der Oberklasse Problemtype -
  die kann ich noch in specific.C unterbringen.

  Hier können beliebige benötigte Subroutinen eingebunden werden -
  damit stören sie nicht im Hauptcode
 */

class LinearWaveDecay: public ProblemType {
public:
	LinearWaveDecay (Data &);
	virtual void init_fields(Data &, int [], int []);
	virtual void writePhystestInfo(Data &);
private:
	double get_FluctData(Data &);
	double get_FluctDataDirect(Data &);
	void set_shearAlfven(Data &, int ix, int iy, int iz);
	void set_fastMode(Data &, int ix, int iy, int iz);
	void set_slowMode(Data &, int ix, int iy, int iz);
	void set_soundMode(Data &, int ix, int iy, int iz);
	void set_ShearFlow(Data &, int ix, int iy, int iz);
	void set_Background(Data &, int ix, int iy, int iz);
	void WriteToH5(Hdf5Stream &h5out);
	NumArray<REAL> fields;
	std::vector<double> timeVelo, rmsVelo;
	REAL k_x, k_y;
	REAL v_amp;
	REAL B0, rho0, p0;
	REAL len0;
	REAL cAlfven, cFast, cSlow, cSound;
	REAL cA_k;
	int wave_type;
};

LinearWaveDecay::LinearWaveDecay(Data &gdata)  : ProblemType(gdata)
{
#if (GEOM == CARTESIAN)
	name = " Alfven wing test ";
#elif (GEOM == CYLINDRICAL) 
	name = " Alfven wing test in cylindrical geometry ";
#elif (GEOM == SPHERICAL)
	name = " Alfven wing test in spherical geometry ";
#endif

	this->v_amp = 0.05;
	this->v_amp = value((char*)"Initial_Amplitude");
	this->wave_type = value((char*)"Wave_Type"); // 1 = shear Alvfen, 2 = fast, 3 = slow

	// set wave number
	this->k_x = 2.*M_PI/gdata.Lx[0];
	this->k_y = 2.*M_PI/gdata.Lx[1];
  
	// Set correct size for array holding disturbance
	fields.resize(8);

	B0 = 1.;
	if(value_exists("BackgroundDensity")) {
		rho0 = value((char*)"BackgroundDensity");
	} else {
		rho0 = 1.;
	}
	p0 = 1.;

	cSound = sqrt(gamma*p0/rho0);
	cAlfven = sqrt(B0*B0/(rho0));

	if(value_exists("Alfven_Speed")) {
		cAlfven = value((char*)"Alfven_Speed");
		B0 = sqr(cAlfven)*rho0;
	}
	if(value_exists("Sound_Speed")) {
		cSound  = value((char*)"Sound_Speed");
		p0 = sqr(cSound)*rho0/gamma;
	}
	len0 = gdata.Lx[0];

	cA_k = cAlfven/sqrt(2.);
	REAL theta = M_PI/4;
	cA_k = cAlfven*cos(theta);
	cFast = sqrt(cSound*cSound + cAlfven*cAlfven);
	cFast = 0.5*(sqr(cSound) + sqr(cAlfven) + sqrt(sqr(sqr(cSound) + sqr(cAlfven)) - 4*sqr(cSound)*sqr(cAlfven)*sqr(cos(theta))));
	cFast = sqrt(cFast);
	cSlow = 0.5*(sqr(cSound) + sqr(cAlfven) - sqrt(sqr(sqr(cSound) + sqr(cAlfven)) - 4*sqr(cSound)*sqr(cAlfven)*sqr(cos(theta))));
	cSlow = sqrt(cSlow);

	if(gdata.rank==0) {
		cout << " c_sound " << cSound << endl;
		cout << " c_Alfven " << cAlfven << endl;
		cout << " c_Alfven in k direction " << cA_k << endl;
		cout << " fast mode " << cFast << endl;
		cout << " slow mode " << cSlow << endl;
	}

//	cF =
//	exit(3);

//	wave_type = 1; // shear Alfven waves
}





void LinearWaveDecay::init_fields(Data &gdata, int ibeg[3], int iend[3])
{

	for(int iz=ibeg[2]; iz<=iend[2]; ++iz) {
		for(int iy=ibeg[1]; iy<=iend[1]; ++iy) {
			for(int ix=ibeg[0]; ix<=iend[0]; ++ix) {

				//					gdata.om[q](ix,iy,iz) = 1.;
				if(wave_type==1) {
					set_shearAlfven(gdata, ix,iy,iz);
				} else if (wave_type==2) {
					set_fastMode(gdata, ix, iy, iz);
				} else if (wave_type==3) {
					set_slowMode(gdata, ix, iy, iz);
				} else if (wave_type==4) {
					set_soundMode(gdata, ix, iy, iz);
				} else if (wave_type==5) {
					set_ShearFlow(gdata, ix, iy, iz);
				}

				gdata.om[q_rho](ix,iy,iz) = rho0 + fields(q_rho);
				gdata.om[q_sx](ix,iy,iz) = fields(q_sx);
				gdata.om[q_sy](ix,iy,iz) = fields(q_sy);
				gdata.om[q_sz](ix,iy,iz) = fields(q_sz);
				gdata.om[q_Bx](ix,iy,iz) = fields(q_Bx);
				gdata.om[q_By](ix,iy,iz) = fields(q_By);
				REAL yPos = gdata.getEdgL_y(iy);
				gdata.om[q_Bz](ix,iy,iz) = B0*yPos + fields(q_Bz);
				gdata.om[q_Eges](ix,iy,iz) = p0;


			}
		}
	}

	if(gdata.rank==0) {
	cout << " A_x ";
	cout << gdata.om[q_Bx](12,-2,0) << " ";
	cout << gdata.om[q_Bx](12,-1,0) << " ";
	cout << gdata.om[q_Bx](12, 0,0) << " ";
	cout << gdata.om[q_Bx](12, 1,0) << " ";
	cout << endl;
	cout << " A_y ";
	cout << gdata.om[q_By](12,-2,0) << " ";
	cout << gdata.om[q_By](12,-1,0) << " ";
	cout << gdata.om[q_By](12, 0,0) << " ";
	cout << gdata.om[q_By](12, 1,0) << " ";
	cout << endl;
	cout << " A_z ";
	cout << gdata.om[q_Bz](12,-3,0) << " ";
	cout << gdata.om[q_Bz](12,-2,0) << " ";
	cout << gdata.om[q_Bz](12,-1,0) << " ";
	cout << gdata.om[q_Bz](12, 0,0) << " ";
	cout << gdata.om[q_Bz](12, 1,0) << " ";
	cout << gdata.om[q_Bz](12, 2,0) << " ";
	cout << endl;
	}
	if(gdata.rank==1) {
		cout << " A_x1 ";
		cout << gdata.om[q_Bx](12,gdata.mx[1]-3,0) << " ";
		cout << gdata.om[q_Bx](12,gdata.mx[1]-2,0) << " ";
		cout << gdata.om[q_Bx](12,gdata.mx[1]-1,0) << " ";
		cout << gdata.om[q_Bx](12,gdata.mx[1]+0,0) << " ";
		cout << gdata.om[q_Bx](12,gdata.mx[1]+1,0) << " ";
		cout << gdata.om[q_Bx](12,gdata.mx[1]+2,0) << " ";
		cout << endl;
		cout << " A_y1 ";
		cout << gdata.om[q_By](12,gdata.mx[1]-3,0) << " ";
		cout << gdata.om[q_By](12,gdata.mx[1]-2,0) << " ";
		cout << gdata.om[q_By](12,gdata.mx[1]-1,0) << " ";
		cout << gdata.om[q_By](12,gdata.mx[1]+0,0) << " ";
		cout << gdata.om[q_By](12,gdata.mx[1]+1,0) << " ";
		cout << gdata.om[q_By](12,gdata.mx[1]+2,0) << " ";
		cout << endl;
		cout << " A_z1 ";
		cout << gdata.om[q_Bz](12,gdata.mx[1]-3,0) << " ";
		cout << gdata.om[q_Bz](12,gdata.mx[1]-2,0) << " ";
		cout << gdata.om[q_Bz](12,gdata.mx[1]-1,0) << " ";
		cout << gdata.om[q_Bz](12,gdata.mx[1]+0,0) << " ";
		cout << gdata.om[q_Bz](12,gdata.mx[1]+1,0) << " ";
		cout << gdata.om[q_Bz](12,gdata.mx[1]+2,0) << " ";
		cout << gdata.om[q_Bz](12,gdata.mx[1]+3,0) << " ";
		cout << endl;
	}

	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");


}


void LinearWaveDecay::set_shearAlfven(Data &gdata, int ix, int iy, int iz) {
	/*!
	 * Shear Alfven wave test taken from
	 * RyuEtAl1995
	 * the setup is given in Eq. (3.6) of that paper
	 * */
	REAL xPos = gdata.getCen_x(ix);
	REAL yPos = gdata.getCen_y(iy);


	REAL del_vz = v_amp*cA_k*sin(k_x*xPos + k_y*yPos);
	REAL del_rho = 0.;
	REAL del_p   = 0.;
	REAL del_vx  = 0.;
	REAL del_vy  = 0.;
	REAL del_Ax  = 0.;
	REAL del_Ay  = 0.;
	REAL del_Az  = 0.;

	fields(q_rho) = del_rho;
	fields(q_sx) = del_vx;
	fields(q_sy) = del_vy;
	fields(q_sz) = del_vz;
	fields(q_Bx) = del_Ax;
	fields(q_By) = del_Ay;
	fields(q_Bz) = del_Az;
	fields(q_Eges) = del_p;

}


void LinearWaveDecay::set_fastMode(Data &gdata, int ix, int iy, int iz) {
	/*!
	 * fast-mode wave test taken from
	 * RyuEtAl1995
	 * the setup is given in Eq. (3.10) of that paper
	 * */
	REAL xPos = gdata.getCen_x(ix);
	REAL yPos = gdata.getCen_y(iy);


	REAL cA = 1.;
	REAL del_vx = v_amp*cFast*sin(k_x*xPos + k_y*yPos);
	REAL del_vy = (sqr(cFast)/(sqr(cFast) - sqr(cAlfven)))*v_amp*cFast*sin(k_x*xPos + k_y*yPos);
	REAL del_rho = 0.;
	REAL del_p   = 0.;
	REAL del_vz  = 0.;
	REAL del_Ax  = 0.;
	REAL del_Ay  = 0.;
	REAL del_Az  = 0.;

	fields(q_rho) = del_rho;
	fields(q_sx) = del_vx;
	fields(q_sy) = del_vy;
	fields(q_sz) = del_vz;
	fields(q_Bx) = del_Ax;
	fields(q_By) = del_Ay;
	fields(q_Bz) = del_Az;
	fields(q_Eges) = del_p;

}
void LinearWaveDecay::set_slowMode(Data &gdata, int ix, int iy, int iz) {
	/*!
	 * slow-mode wave test taken from
	 * RyuEtAl1995
	 * the setup is given in Eq. (3.11) of that paper
	 * */
	REAL xPos = gdata.getCen_x(ix);
	REAL yPos = gdata.getCen_y(iy);


	REAL cA = 1.;
	REAL del_vx = v_amp*cSlow*sin(k_x*xPos + k_y*yPos);
	REAL del_vy = (sqr(cSlow)/(sqr(cSlow) - sqr(cAlfven)))*v_amp*cSlow*sin(k_x*xPos + k_y*yPos);
	REAL del_rho = 0.;
	REAL del_p   = 0.;
	REAL del_vz  = 0.;
	REAL del_Ax  = 0.;
	REAL del_Ay  = 0.;
	REAL del_Az  = 0.;

	fields(q_rho) = del_rho;
	fields(q_sx) = del_vx;
	fields(q_sy) = del_vy;
	fields(q_sz) = del_vz;
	fields(q_Bx) = del_Ax;
	fields(q_By) = del_Ay;
	fields(q_Bz) = del_Az;
	fields(q_Eges) = del_p;
}



void LinearWaveDecay::set_soundMode(Data &gdata, int ix, int iy, int iz) {
	/*! Test setup according to:
	 * RyuGoodman1994ApJ422_269
	 * this setup represents a decaying sound wave given in Eq (51) in that paper
	 * */
	REAL xPos = gdata.getCen_x(ix);
	REAL yPos = gdata.getCen_y(iy);

	REAL del_vx = v_amp*cSound*sin(k_x*xPos + k_y*yPos);
	REAL del_vy = del_vx;
	REAL del_rho = 0.;
	REAL del_p   = 0.;
	REAL del_vz  = 0.;
	REAL del_Ax  = 0.;
	REAL del_Ay  = 0.;
	REAL del_Az  = 0.;

	fields(q_rho) = del_rho;
	fields(q_sx) = del_vx;
	fields(q_sy) = del_vy;
	fields(q_sz) = del_vz;
	fields(q_Bx) = del_Ax;
	fields(q_By) = del_Ay;
	fields(q_Bz) = del_Az;
	fields(q_Eges) = del_p;
}



void LinearWaveDecay::set_ShearFlow(Data &gdata, int ix, int iy, int iz) {
	/*! Test setup according to:
	 * RyuGoodman1994ApJ422_269
	 * this setup represents a planar shear flow, where the gradual decay can be investigated
	 * the specific setup is given in Eq. (54) in that paper
	 * */


	REAL xPos = gdata.getCen_x(ix);
	REAL yPos = gdata.getCen_y(iy);

	REAL del_vx = v_amp*cSound*sin(k_x*(xPos - yPos));
	REAL del_vy = del_vx;
	REAL del_rho = 0.;
	REAL del_p   = 0.;
	REAL del_vz  = 0.;
	REAL del_Ax  = 0.;
	REAL del_Ay  = 0.;
	REAL del_Az  = 0.;

	fields(q_rho) = del_rho;
	fields(q_sx) = del_vx;
	fields(q_sy) = del_vy;
	fields(q_sz) = del_vz;
	fields(q_Bx) = del_Ax;
	fields(q_By) = del_Ay;
	fields(q_Bz) = del_Az;
	fields(q_Eges) = del_p;

}



double LinearWaveDecay::get_FluctData(Data &gdata) {
	double rms_vx = 0.;
	double rms_vy = 0.;
	double rms_vz = 0.;
	double num = 0.;
	for( int ix=0; ix<=gdata.mx[0]; ++ix) {
		for( int iy=0; iy<=gdata.mx[1]; ++iy) {
			num += 1.;
			rms_vx+= sqr(gdata.om[q_sx](ix,iy,0));
			rms_vy+= sqr(gdata.om[q_sy](ix,iy,0));
			rms_vz+= sqr(gdata.om[q_sz](ix,iy,0));
		}
	}
#ifdef parallel
	double numGlobal, rms_vxGlobal, rms_vyGlobal, rms_vzGlobal;
	MPI_Allreduce(&num, &numGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Allreduce(&rms_vx, &rms_vxGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Allreduce(&rms_vy, &rms_vyGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Allreduce(&rms_vz, &rms_vzGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Barrier(gdata.comm3d);
	num = numGlobal;
	rms_vx = rms_vxGlobal;
	rms_vy = rms_vyGlobal;
	rms_vz = rms_vzGlobal;
#endif

	rms_vx /= num;
//	rms_vx = sqrt(rms_vx);
	rms_vy /= num;
//	rms_vy = sqrt(rms_vy);
	rms_vz /= num;
//	rms_vz = sqrt(rms_vz);

	double rms_velo = 0.;
	if(wave_type==1) {
		rms_velo = rms_vz;
	} else if (wave_type==2) {
		rms_velo = sqr(cFast)*rms_vx + sqr(sqr(cFast)/(sqr(cFast) - sqr(cAlfven)))*rms_vy;
	} else if (wave_type==3) {
		rms_velo = sqr(cSlow)*rms_vx + sqr(sqr(cSlow)/(sqr(cSlow) - sqr(cAlfven)))*rms_vy;
	} else if (wave_type==4) {
		rms_velo = rms_vx;
	} else if (wave_type==5) {
		rms_velo = rms_vx;
	} else {
		cerr << " Error: no such wave type implemented - Exiting " << endl;
		exit(3);
	}
	return rms_velo;
}


double LinearWaveDecay::get_FluctDataDirect(Data &gdata) {
	double amp_vx(0.), amp_vy(0.), amp_vz(0.);
	double amplitude(0.);
	double num(0.);

	for( int ix=0; ix<=gdata.mx[0]; ++ix) {
		for( int iy=0; iy<=gdata.mx[1]; ++iy) {
			if(wave_type==1) {
				set_shearAlfven(gdata, ix,iy,0);
				double val_init = fields(q_sz);
				if(std::abs(val_init)>1.e-12) {
					amp_vz += sqr(gdata.om[q_sz](ix,iy,0)/val_init);
					num += 1.;
				}
			} else {
				if (wave_type==2) {
					set_fastMode(gdata, ix, iy, 0);
				} else if (wave_type==3) {
					set_slowMode(gdata, ix, iy, 0);
				} else if (wave_type==4){
					set_soundMode(gdata, ix, iy, 0);
				} else {
					set_ShearFlow(gdata, ix, iy, 0);
				}
				double val_init = fields(q_sy);
				if(std::abs(val_init)>1.e-12) {
					amp_vy += sqr(gdata.om[q_sy](ix,iy,0)/val_init);
					num += 1.;
				}
			}
//			if(iy==12) {
//				cout << " vals " << amp_vy << " " << fields(q_sy) << " " << gdata.om[q_sy](ix,iy,0) << " " << num << endl;
//			}
		}
	}

#ifdef parallel
	double numGlobal, amp_vxGlobal, amp_vyGlobal, amp_vzGlobal;
	MPI_Allreduce(&num, &numGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Allreduce(&amp_vx, &amp_vxGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Allreduce(&amp_vy, &amp_vyGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Allreduce(&amp_vz, &amp_vzGlobal, 1, MPI_DOUBLE, MPI_SUM, gdata.comm3d);
	MPI_Barrier(gdata.comm3d);
	num = numGlobal;
	amp_vx = amp_vxGlobal;
	amp_vy = amp_vyGlobal;
	amp_vz = amp_vzGlobal;
#endif
	amp_vx /= num;
	amp_vy /= num;
	amp_vz /= num;

	if(wave_type==1) {
		amplitude = amp_vz;
	} else {
		amplitude = amp_vy;
	}

//	cout << " amp: " << amplitude << endl;
//	exit(3);

	return amplitude;

}

void LinearWaveDecay::writePhystestInfo(Data &gdata) {
	// Here we use this to compute the fluctuation statistics at each relevant step


	double rms_velo = get_FluctData(gdata);
//	double rms_velo = get_FluctDataDirect(gdata);

	rmsVelo.push_back(sqrt(rms_velo));
	timeVelo.push_back(gdata.time);

}



void LinearWaveDecay::WriteToH5(Hdf5Stream &h5out) {

	  // Indicate that we are using a test case
	  h5out.AddGlobalAttr("TestCase", "Linear wave decay");

	  // supply normalised values necessary for Parker solution
	  // Use own group for that purpose
	  hid_t group = h5out.AddGroup("Data/WavePars");

	  // supply wave type
	  h5out.AddGlobalAttr("Wave_Type", wave_type, group);

	  // suppy initial amplitude
	  h5out.AddGlobalAttr("Initial_Amplitude", v_amp, group);

	  // supply normalised speed of sound
	  h5out.AddGlobalAttr("cs", cSound, group);

	  // supply normalised Alfven speed
	  h5out.AddGlobalAttr("cA", cAlfven, group);

	  // supply Alfven wave speed
	  h5out.AddGlobalAttr("cWave_Alfven", cA_k, group);

	  // supply fast mode speed
	  h5out.AddGlobalAttr("cWave_Fast", cFast, group);

	  // supply slow mode speed
	  h5out.AddGlobalAttr("cWave_Slow", cSlow, group);

	  // Supply additional extent of domain
	  h5out.AddGlobalAttr("Len", len0, group);

	  // Background density
	  h5out.AddGlobalAttr("BackgroundDensity", rho0, group);


	  // Add full array of times and wave amplitudes
	  NumArray<float> outData;
	  outData.resize(rmsVelo.size());
	  for(int i_step=0; i_step<rmsVelo.size(); ++i_step) {
		  outData(i_step) = timeVelo[i_step];
	  }
	  h5out.WriteNumArray("times",outData,group);
	  for(int i_step=0; i_step<rmsVelo.size(); ++i_step) {
		  outData(i_step) = rmsVelo[i_step];
	  }
	  h5out.WriteNumArray("rmsVelo",outData,group);

	  // Now close the new group
	  h5out.CloseGroup(group);

	}



#endif
