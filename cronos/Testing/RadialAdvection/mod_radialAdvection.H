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

class RadialAdvectionTest: public ProblemType {
public:
	RadialAdvectionTest (Data &);
	virtual void init_fields(Data &, int [], int []);
	virtual void src_User(Data &, NumMatrix<REAL,3> [], NumMatrix<REAL,3> []);
//	virtual void writePhystestInfo(Data &);
	virtual void bc_User(Data &, NumMatrix<REAL,3> &, int, int, int, int);
  virtual void computePhystest(Data &);
private:
	void WriteToH5(Hdf5Stream &h5out);
	double get_rmsError(Data &);
	double get_velo(Data &, int ir);
	double get_velo_x(Data &, int ix);
	double get_velo_y(Data &, int iy);
	double v0;
};

RadialAdvectionTest::RadialAdvectionTest(Data &gdata)  : ProblemType(gdata)
{
#if (GEOM == CYLINDRICAL)
	name = " radial advection test ";
#elif (GEOM == CARTESIAN)
	name = " radial advection test on Cartisian grid ";
#else
	cerr << " Not implemented, yet " << endl;
	exit(3);
#endif

	this->v0 = value((char*)"Velocity_scale");
  
}





void RadialAdvectionTest::init_fields(Data &gdata, int ibeg[3], int iend[3])
{

	for(int iz=ibeg[2]; iz<=iend[2]; ++iz) {
		for(int iy=ibeg[1]; iy<=iend[1]; ++iy) {
			for(int ix=ibeg[0]; ix<=iend[0]; ++ix) {
			  // The velocity is linear in all cases ->
			  // should not determine the order
				gdata.om[q_rho](ix,iy,iz) = rho0;
#if GEOM == CYLINDRICAL
				gdata.om[q_sx](ix,iy,iz) = get_velo(gdata, ix);
				gdata.om[q_sy](ix,iy,iz) = 0.;
				gdata.om[q_sz](ix,iy,iz) = 0.;
#else
				gdata.om[q_sx](ix,iy,iz) = get_velo_x(gdata, ix);
				gdata.om[q_sy](ix,iy,iz) = get_velo_y(gdata, iy);
				gdata.om[q_sz](ix,iy,iz) = 0.;
#endif

			}
		}
	}


}

double RadialAdvectionTest::get_velo(Data &gdata, int ir) {
	double r_cyl = gdata.getCen_x(ir);
	return v0*r_cyl;

}


double RadialAdvectionTest::get_velo_x(Data &gdata, int ix) {
	double xPos = gdata.getCen_x(ix);
	double v_x = v0*xPos;
	return v_x;

}

double RadialAdvectionTest::get_velo_y(Data &gdata, int iy) {
	double yPos = gdata.getCen_y(iy);
	double v_y = v0*yPos;
	return v_y;

}


//void RadialAdvectionTest::writePhystestInfo(Data &gdata) {
//	// Here we use this to compute the fluctuation statistics at each relevant step
//
//
//	double rms_velo = get_FluctData(gdata);
////	double rms_velo = get_FluctDataDirect(gdata);
//
//	rmsVelo.push_back(sqrt(rms_velo));
//	timeVelo.push_back(gdata.time);
//
//}



void RadialAdvectionTest::src_User(Data &gdata, NumMatrix<REAL,3> nom[N_OMINT],
                         NumMatrix<REAL,3> nom_user[N_OMINT_USER])
{
  for(int iz=0; iz<=gdata.mx[2]; ++iz) {
    for(int iy=0; iy<=gdata.mx[1]; ++iy) {
      for(int ix=0; ix<=gdata.mx[0]; ++ix) {

    	  // Using pure advection setup
    	  nom[q_sx](ix,iy,iz) = 0.;
#if GEOM == CYLINDRICAL
    	  nom[q_sx](ix,iy,iz) = nom[q_rho](ix,iy,iz)*get_velo(gdata,ix);
    	  nom[q_sy](ix,iy,iz) = 0.;
    	  nom[q_sz](ix,iy,iz) = 0.;
#else
    	  nom[q_sx](ix,iy,iz) = nom[q_rho](ix,iy,iz)*get_velo_x(gdata,ix);
    	  nom[q_sy](ix,iy,iz) = nom[q_rho](ix,iy,iz)*get_velo_y(gdata,iy);
    	  nom[q_sz](ix,iy,iz) = 0.;
#endif
//    	  gdata.om[q_sx](ix,iy,iz) = get_velo(gdata,ix);

      }
    }
  }
  //  get_rmsError(gdata);
}


void RadialAdvectionTest::bc_User(Data &gdata, NumMatrix<REAL,3> &omb,
		int dir, int top, int q, int rim)
{
	/*
      Include fixed boundary conditions for the z-Axis:
  Ghost cells will be filled with the initial values
	 */

#if(GEOM==CYLINDRICAL)
	if(dir == 0) {

		if (top == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; ++k) {
				for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
					for (int i = -rim; i <0; ++i) {
						if(q==1) {
							omb(i,j,k) = get_velo(gdata,i);
						} else if (q==0) {
							omb(i,j,k) = omb(5,j,k);
						} else {
							omb(i,j,k) = 0.;
						}
					}
				}
			}

		} else {
			for (int k = -rim; k <= gdata.mx[2]+rim; ++k) {
				for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
					for (int i = gdata.mx[0]+1; i <=gdata.mx[0]+rim; ++i) {
						if(q==1) {
//							omb(i,j,k) = get_velo(gdata,i);
							omb(i,j,k) = 2*omb(i-1,j,k) - omb(i-2,j,k);
							// if(j==0 && k==0) {
							//   cout << " Extrapol " << i << " ";
							//   cout << omb(i-1,j,k) << " " << omb(i-2,j,k) << " ";
							//   cout << omb(i,j,k) << endl;
							// }
						} else if (q==0) {
							omb(i,j,k) = omb(gdata.mx[0]-3,j,k);
						} else {
							omb(i,j,k) = 0.;
						}
					}
				}
			}

		}

	}
#else
	if(dir == 0) {

		if (top == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; ++k) {
				for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
					for (int i = -rim; i <0; ++i) {
						if(q==1) {
							omb(i,j,k) = get_velo_x(gdata,i);
						} else if (q==2) {
							omb(i,j,k) = get_velo_y(gdata,j);
						} else if (q==0) {
							omb(i,j,k) = omb(5,j,k);
						} else {
							omb(i,j,k) = 0.;
						}
					}
				}
			}

		} else {
			for (int k = -rim; k <= gdata.mx[2]+rim; ++k) {
				for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
					for (int i = gdata.mx[0]+1; i <=gdata.mx[0]+rim; ++i) {
						if(q==1) {
//							omb(i,j,k) = get_velo(gdata,i);
							omb(i,j,k) = 2*omb(i-1,j,k) - omb(i-2,j,k);
							// if(j==0 && k==0) {
							//   cout << " Extrapol " << i << " ";
							//   cout << omb(i-1,j,k) << " " << omb(i-2,j,k) << " ";
							//   cout << omb(i,j,k) << endl;
							// }
						} else if (q==2) {
							omb(i,j,k) = get_velo_y(gdata,j);
						} else if (q==0) {
							omb(i,j,k) = omb(gdata.mx[0]-3,j,k);
						} else {
							omb(i,j,k) = 0.;
						}
					}
				}
			}

		}

	} else if (dir==1) {
		if (top == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; ++k) {
				for (int j = -rim; j <0; ++j) {
					for (int i = -rim; i <= gdata.mx[0]+rim; i++) {
						if(q==1) {
							omb(i,j,k) = get_velo_x(gdata,i);
						} else if (q==2) {
							omb(i,j,k) = get_velo_y(gdata,j);
						} else if (q==0) {
							omb(i,j,k) = omb(5,j,k);
						} else {
							omb(i,j,k) = 0.;
						}
					}
				}
			}

		} else {
			for (int k = -rim; k <= gdata.mx[2]+rim; ++k) {
				for (int j = gdata.mx[1]+1; j <=gdata.mx[1]+rim; ++j) {
					for (int i = -rim; i <= gdata.mx[0]+rim; i++) {
						if(q==1) {
//							omb(i,j,k) = get_velo(gdata,i);
							omb(i,j,k) = 2*omb(i-1,j,k) - omb(i-2,j,k);
							// if(j==0 && k==0) {
							//   cout << " Extrapol " << i << " ";
							//   cout << omb(i-1,j,k) << " " << omb(i-2,j,k) << " ";
							//   cout << omb(i,j,k) << endl;
							// }
						} else if (q==2) {
							omb(i,j,k) = get_velo_y(gdata,j);
						} else if (q==0) {
							omb(i,j,k) = omb(gdata.mx[0]-3,j,k);
						} else {
							omb(i,j,k) = 0.;
						}
					}
				}
			}

		}

	}
#endif
}


void RadialAdvectionTest::computePhystest(Data &gdata) {
  // Get rms error:
  get_rmsError(gdata);
}


double RadialAdvectionTest::get_rmsError(Data &gdata) {
	double ana = rho0*exp(-2.*v0*gdata.time);
	double err_rms = 0.;
	double ave_val = 0.;
	double num = 0.;

	for(int iz=0; iz<=gdata.mx[2]; ++iz) {
		for(int iy=0; iy<=gdata.mx[1]; ++iy) {
			for(int ix=0; ix<=gdata.mx[0]; ++ix) {
				double diff = ana - gdata.om[q_rho](ix,iy,iz);
				err_rms += sqr(diff);
				ave_val += gdata.om[q_rho](ix,iy,iz);
				num++;
			}
		}
	}



	ave_val /= num;
	err_rms /= num;
	err_rms = sqrt(err_rms);
	if(gdata.rank==0) {
		cout << " rms deviation: " << err_rms << "-> relative " << err_rms/ave_val << endl;
	}

	num = 0.;
	double standDev=0;
	// Also compute related statistics
	for(int iz=0; iz<=gdata.mx[2]; ++iz) {
		for(int iy=0; iy<=gdata.mx[1]; ++iy) {
			for(int ix=0; ix<=gdata.mx[0]; ++ix) {
				double diff = gdata.om[q_rho](ix,iy,iz) - ave_val;
				standDev += sqr(diff);
				num++;
			}
		}
	}
	standDev = sqrt(standDev)/(num-1);

	if(gdata.rank==0) {
		cout << " rms value: " << ave_val << " +/- " << standDev << endl;
	}


	return err_rms;
}


void RadialAdvectionTest::WriteToH5(Hdf5Stream &h5out) {

	  // Indicate that we are using a test case
	  h5out.AddGlobalAttr("TestCase", "Radial advection");

	  // supply normalised values necessary for Parker solution
	  // Use own group for that purpose
	  hid_t group = h5out.AddGroup("Data/WavePars");

	  // supply wave type
	  h5out.AddGlobalAttr("Velocity scale", v0, group);

	  // Now close the new group
	  h5out.CloseGroup(group);

	}



#endif
