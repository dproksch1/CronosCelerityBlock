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

class AlfvenWingsSpherical: public ProblemType {
public:
	AlfvenWingsSpherical(Data &);
	virtual void init_fields(Data &, int [], int []);
	virtual void bc_User(Data &, NumMatrix<REAL,3> &, int, int, int, int);
	virtual void src_User(Data &, NumMatrix<REAL,3> []);
private:
	REAL GetConstFields(Data &, int q, int ii, int jj, int kk);
	REAL Get_x1Vec(Data &, REAL valx, REAL valy, REAL valz,
	               int ii, int jj, int kk);
	REAL Get_x2Vec(Data &, REAL valx, REAL valy, REAL valz,
	               int ii, int jj, int kk);
	REAL Get_x3Vec(Data &, REAL valx, REAL valy, REAL valz,
	               int ii, int jj, int kk);

	REAL Get_x1Vec(Data &, REAL valx, REAL valy, REAL valz, REAL ipos[3]);
	REAL Get_x2Vec(Data &, REAL valx, REAL valy, REAL valz, REAL ipos[3]);
	REAL Get_x3Vec(Data &, REAL valx, REAL valy, REAL valz, REAL ipos[3]);


	REAL Get_x1Vec(Data &, REAL valx, REAL valy, REAL valz,
	               REAL ii, REAL jj, REAL kk);
	REAL Get_x2Vec(Data &, REAL valx, REAL valy, REAL valz,
	               REAL ii, REAL jj, REAL kk);
	REAL Get_x3Vec(Data &, REAL valx, REAL valy, REAL valz,
	               REAL ii, REAL jj, REAL kk);

	REAL Ax(Data &, REAL ipos[3]);
	REAL Ay(Data &, REAL ipos[3]);
	REAL Az(Data &, REAL ipos[3]);

	REAL Ax(Data &, REAL ii, REAL jj, REAL kk);
	REAL Ay(Data &, REAL ii, REAL jj, REAL kk);
	REAL Az(Data &, REAL ii, REAL jj, REAL kk);

	REAL Ax(Data &, int ii, int jj, int kk);
	REAL Ay(Data &, int ii, int jj, int kk);
	REAL Az(Data &, int ii, int jj, int kk);

	REAL B0(Data &, int ii, int jj, int kk);
	REAL B1(Data &, int ii, int jj, int kk);
	REAL B2(Data &, int ii, int jj, int kk);
	REAL vBack[3], BBack[3];
	REAL xpos[3];
	REAL eth0;
};


AlfvenWingsSpherical::AlfvenWingsSpherical(Data &gdata)
{
#if (GEOM == CARTESIAN)
	name = " Alfven wing test ";
#elif (GEOM == CYLINDRICAL) 
	name = " Alfven wing test in cylindrical geometry ";
#elif (GEOM == SPHERICAL)
	name = " Alfven wing test in spherical geometry ";
#endif

#if (GEOM != 3)
	cerr << " Error: Supposed to be a spherical problem ";
	cerr << endl;
	exit(2);
#endif

	// Standard configuration
	this->vBack[0] = 0.;
	this->vBack[1] = 0.;
	this->vBack[2] = 1.;

	this->BBack[0] = 1.;
	this->BBack[1] = 0.;
	this->BBack[2] = 0.;

	this->xpos[0]=value((char*)"xPosition");
	this->xpos[1]=value((char*)"yPosition");
	this->xpos[2]=value((char*)"zPosition");

	this->eth0 = cs2*rho0/(gamma-1.);
  
}



void AlfvenWingsSpherical::src_User(Data &gdata,
                                    NumMatrix<REAL,3> nom[N_OMINT]) {

	// 	REAL rhomin = value((char*)"rhomin");
	REAL timeCorr(1.e3);
// 	REAL tau = timeCorr*(std::min(10.*gdata.time,1.));
	REAL tau = std::min(10.*gdata.time,1.);

	for(int k=0; k<=gdata.mx[2]; ++k) {
		for(int j=0; j<=gdata.mx[1]; ++j) {
			for(int i=0; i<=gdata.mx[0]; ++i) {

				REAL xx = gdata.get_xCart(i,j,k);
				REAL yy = gdata.get_yCart(i,j,k);
				REAL zz = gdata.get_zCart(i,j,k);

				REAL xdist = (xx - xpos[0]);
				REAL ydist = (yy - xpos[1]);
				REAL zdist = (zz - xpos[2]);

				// 	      om[0](i,j,k) = std::max(om[0](i,j,k),rhomin);
				double dist = sqrt(sqr(xdist) + sqr(ydist) + sqr(zdist));
				double scaling = 1-tau*(1-tanh(4.*std::max(4.*dist-1.,0.)));

// 				REAL change = tau*(1-tanh(4.*std::max(2.*dist-1.,0.)));
// 				for (int q=1; q<=3; ++q) {
// 					nom[q](i,j,k) = change*gdata.om[q](i,j,k);
// 				}

				// Scaling of velocity:
				for (int q = 1; q <= 3; q++) {
					gdata.om[q](i,j,k) *= scaling;
				}

			}
		}
	}




}



void AlfvenWingsSpherical::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{

	for(int q=0; q<8; ++ q) {
		for(int k=ibeg[2]; k<=iend[2]; ++k) {
			for(int j=ibeg[1]; j<=iend[1]; ++j) {
				for(int i=ibeg[0]; i<=iend[0]; ++i) {

					gdata.om[q](i,j,k) = GetConstFields(gdata, q,i,j,k);
				  
				}
			}
		}
	}

	int i=10;
	int j=11;
	int k=12;

	REAL f_geom = 1./(gdata.h0(i,j,k)*gdata.h1(i,j,k)*gdata.h2(i,j,k));
  
	double dxBx = (gdata.h1(i+0.5,j,k)*gdata.h2(i+0.5,j,k)*gdata.om[4](i,j,k) -
	               gdata.h1(i-0.5,j,k)*gdata.h2(i-0.5,j,k)*gdata.om[4](i-1,j,k))/gdata.dx[0]*f_geom;
  
	double dyBy = (gdata.h0(i,j+0.5,k)*gdata.h2(i,j+0.5,k)*gdata.om[5](i,j,k) - 
	               gdata.h0(i,j-0.5,k)*gdata.h2(i,j-0.5,k)*gdata.om[5](i,j-1,k))/gdata.dx[1]*f_geom;

	double dzBz = (gdata.h0(i,j,k+0.5)*gdata.h1(i,j,k+0.5)*gdata.om[6](i,j,k) -
	               gdata.h0(i,j,k-0.5)*gdata.h1(i,j,k-0.5)*gdata.om[6](i,j,k-1))/gdata.dx[2]*f_geom;

	//   cout << " Test: " << dxBx + dyBy + dzBz << endl;
	//   cout << "       ";
	//   cout << dxBx << " ";
	//   cout << dyBy << " ";
	//   cout << dzBz << " ";
	//   cout << endl;

	//   cout << "       ";
	//   cout << gdata.om[4](i,j,k) << " " << gdata.om[4](i-1,j,k) << " ";
	//   cout << gdata.om[5](i,j,k) << " " << gdata.om[5](i,j-1,k) << " ";
	//   cout << gdata.om[6](i,j,k) << " " << gdata.om[6](i,j,k-1) << " ";
	//   cout << endl;
	//   exit(2);
  

	gdata.om[4].rename("B_x");
	gdata.om[5].rename("B_y");
	gdata.om[6].rename("B_z");


}


REAL AlfvenWingsSpherical::GetConstFields(Data &gdata, int q,
                                          int ii, int jj, int kk) {
	if(q==0) {
		// Density
		return rho0;
	} else if (q==1) {
		// v_x / v_r
		return Get_x1Vec(gdata, vBack[0], vBack[1], vBack[2], ii, jj, kk);
	} else if (q==2) {
		// v_y / v_phi / v_theta
		return Get_x2Vec(gdata, vBack[0], vBack[1], vBack[2], ii, jj, kk);
	} else if (q==3) {
		// v_z / v_phi
		return Get_x3Vec(gdata, vBack[0], vBack[1], vBack[2], ii, jj, kk);
	} else if (q==4) {
		// B_x / B_r
		// 		return Get_x1Vec(gdata, BBack[0], BBack[1], BBack[2], ii, jj, kk);
		return B0(gdata, ii, jj, kk);
	} else if (q==5) {
		// B_y / B_theta / B_phi
		// 		return Get_x2Vec(gdata, BBack[0], BBack[1], BBack[2], ii, jj, kk);
		return B1(gdata, ii, jj, kk);
	} else if (q==6) {
		// B_z / B_phi
		// 		return Get_x3Vec(gdata, BBack[0], BBack[1], BBack[2], ii, jj, kk);
		return B2(gdata, ii, jj, kk);
	} else if (q==7) {
		// Thermal energy
		return eth0;
	} else {
		cerr << " No such variable in this case " << endl;
		exit(2);
		return -99;
	}
}


REAL AlfvenWingsSpherical::Get_x1Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     int ii, int jj, int kk) {
	return Get_x1Vec(gdata, valx, valy, valz, static_cast<REAL>(ii),
	                 static_cast<REAL>(jj), static_cast<REAL>(kk));
}


REAL AlfvenWingsSpherical::Get_x1Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     REAL ipos[3]) {
	return Get_x1Vec(gdata, valx, valy, valz, ipos[0], ipos[1], ipos[2]);
}


REAL AlfvenWingsSpherical::Get_x1Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     REAL ii, REAL jj, REAL kk) {
#if (GEOM == CARTESIAN)

	return valx;

#elif (GEOM == CYLINDRICAL)

	REAL phi = gdata.get_y(jj);
	return (valx*cos(phi) - valy*sin(phi));

#elif (GEOM == SPHERICAL)

	REAL theta = gdata.get_y(jj);
	REAL phi  = gdata.get_z(kk);
	return (valx*sin(theta)*cos(phi) +
	        valy*sin(theta)*sin(phi) +
	        valz*cos(theta));

#endif
}


REAL AlfvenWingsSpherical::Get_x2Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     int ii, int jj, int kk) {
	return Get_x2Vec(gdata, valx, valy, valz, static_cast<REAL>(ii),
	                 static_cast<REAL>(jj), static_cast<REAL>(kk));
}

REAL AlfvenWingsSpherical::Get_x2Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     REAL ipos[3]) {
	return Get_x2Vec(gdata, valx, valy, valz, ipos[0], ipos[1], ipos[2]);
}

REAL AlfvenWingsSpherical::Get_x2Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     REAL ii, REAL jj, REAL kk) {

#if (GEOM == CARTESIAN)

	return valy;

#elif (GEOM == CYLINDRICAL)

	REAL phi = gdata.get_y(jj);
	return (valx*sin(phi) + valy*cos(phi));

#elif (GEOM == SPHERICAL)

	REAL theta = gdata.get_y(jj);
	REAL phi   = gdata.get_z(kk);
	return (valx*cos(theta)*cos(phi) +
	        valy*cos(theta)*sin(phi) -
	        valz*sin(theta));

#endif
}


REAL AlfvenWingsSpherical::Get_x3Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     int ii, int jj, int kk) {
	return Get_x3Vec(gdata, valx, valy, valz, static_cast<REAL>(ii),
	                 static_cast<REAL>(jj), static_cast<REAL>(kk));
}

REAL AlfvenWingsSpherical::Get_x3Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     REAL ipos[3]) {
	return Get_x3Vec(gdata, valx, valy, valz, ipos[0], ipos[1], ipos[2]);
}


REAL AlfvenWingsSpherical::Get_x3Vec(Data &gdata,
                                     REAL valx, REAL valy, REAL valz,
                                     REAL ii, REAL jj, REAL kk) {
#if (GEOM == CARTESIAN)

	return valz;

#elif (GEOM == CYLINDRICAL)

	return valz;

#elif (GEOM == SPHERICAL)

	REAL phi = gdata.get_z(kk);
	return -valx*sin(phi) + valy*cos(phi);

#endif
}


// Cartesian Vector-Potential components:

REAL AlfvenWingsSpherical::Ax(Data &gdata, REAL ii, REAL jj, REAL kk) {
	// A_x = By0*zz;
	return BBack[1]*gdata.get_zCart(ii, jj, kk);
}

REAL AlfvenWingsSpherical::Ax(Data &gdata, REAL ipos[3]) {
	return Ax(gdata, ipos[0], ipos[1], ipos[2]);
}

REAL AlfvenWingsSpherical::Ax(Data &gdata, int ii, int jj, int kk) {
	return Ax(gdata, static_cast<REAL>(ii), static_cast<REAL>(jj),
	          static_cast<REAL>(jj));
}


REAL AlfvenWingsSpherical::Ay(Data &gdata, REAL ii, REAL jj, REAL kk) {
	// A_y = Bz0*xx;
	return BBack[2]*gdata.get_xCart(ii, jj, kk);
}

REAL AlfvenWingsSpherical::Ay(Data &gdata, REAL ipos[3]) {
	return Ay(gdata, ipos[0], ipos[1], ipos[2]);
}

REAL AlfvenWingsSpherical::Ay(Data &gdata, int ii, int jj, int kk) {
	return Ay(gdata, static_cast<REAL>(ii), static_cast<REAL>(jj),
	          static_cast<REAL>(jj));
}


REAL AlfvenWingsSpherical::Az(Data &gdata, REAL ii, REAL jj, REAL kk) {
	// A_z = Bx0*yy;
	return BBack[0]*gdata.get_yCart(ii, jj, kk);
}

REAL AlfvenWingsSpherical::Az(Data &gdata, REAL ipos[3]) {
	return Az(gdata, ipos[0], ipos[1], ipos[2]);
}

REAL AlfvenWingsSpherical::Az(Data &gdata, int ii, int jj, int kk) {
	return Az(gdata, static_cast<REAL>(ii), static_cast<REAL>(jj),
	          static_cast<REAL>(jj));
}


// Compute Mag-Field components from vector-Potential (here we use
// corresponding geometry):

 
REAL AlfvenWingsSpherical::B0(Data &gdata, int ii, int jj, int kk) {

#if (GEOM == CARTESIAN)

	REAL pos[3] = {ii+0.5,jj+0.5,1.*kk};
	REAL Bx  = (Az(gdata, pos) - Az(gdata, pos))*gdata.idx[1];

	pos[0] = ii+0.5; pos[1] = 1.*jj; pos[2] = kk+0.5;
	Bx -= (Ay(gdata, pos) - Ay(gdata, pos))*gdata.idx[2];

	return Bx;

#else 

	REAL pos[3] = {ii+0.5,jj+0.5,1.*kk};
	REAL A_3JP = Get_x3Vec(gdata, Ax(gdata, pos), Ay(gdata, pos),
	                       Az(gdata, pos), pos);

	pos[0] = ii+0.5; pos[1] = jj-0.5; pos[2] = 1.*kk;
	REAL A_3JM = Get_x3Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);

	pos[0] = ii+0.5; pos[1] = 1.*jj; pos[2] = kk+0.5;
	REAL A_2KP = Get_x2Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);
	pos[0] = ii+0.5; pos[1] = 1.*jj; pos[2] = kk-0.5;
	REAL A_2KM = Get_x2Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);

#if (GEOM == CYLINDRICAL)

	REAL r_cyl = gdata.get_x(ii+0.5);

	REAL Br = ((A_3JP - A_3JM)*gdata.idx[1]/r_cyl -
	           (A_2KP - A_2KM)*gdata.idx[2]);

#elif (GEOM == SPHERICAL) 

	REAL r_sph = gdata.get_x(ii+0.5);
	REAL theta = gdata.get_y( 1.*jj);

	REAL sintJP = sin(gdata.get_y(jj+0.5));
	REAL sintJM = sin(gdata.get_y(jj-0.5));

	REAL Br = ((A_3JP*sintJP - A_3JM*sintJM)*gdata.idx[1] -
	           (A_2KP        - A_2KM       )*gdata.idx[2])/(r_sph*sin(theta));
	
#endif
	return Br;
#endif	
}


REAL AlfvenWingsSpherical::B1(Data &gdata, int ii, int jj, int kk) {

#if (GEOM == CARTESIAN)

// 	REAL By = ((Ax(1.*ii,jj+0.5,kk+0.5) - Ax(1.*ii,jj+0.5,kk-0.5))*gdata.idx[1] -
// 	           (Az(ii+0.5,jj+0.5,1.*kk) - Az(ii-0.5,jj+0.5,1.*kk))*gdata.idx[2]);
	REAL pos[3] = {1.*ii,jj+0.5,kk+0.5};
	REAL By = (Ax(gdata, pos) - Ax(gdata, pos))*gdata.idx[1];

	pos[0] = ii+0.5; pos[1] = jj+0.5; pos[2] = 1.*kk;
	By -= (Az(gdata, pos) - Az(gdata, pos))*gdata.idx[2];
	return By;

#else 

	REAL pos[3] = {1.*ii,jj+0.5,kk+0.5};
	REAL A_1KP = Get_x1Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);

	pos[0] = 1.*ii; pos[1] = jj+0.5; pos[2] = kk-0.5;
	REAL A_1KM = Get_x1Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);

	pos[0] = ii+0.5; pos[1] = jj+0.5; pos[2] = 1.*kk;
	REAL A_3IP = Get_x3Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);

	pos[0] = ii-0.5; pos[1] = jj+0.5; pos[2] = 1.*kk;
	REAL A_3IM = Get_x3Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);


#if (GEOM == CYLINDRICAL)

	REAL r_cyl = gdata.get_x(ii);

	REAL Bphi = ((A_1KP - A_1KM)*gdata.idx[2] -
	             (A_3IP - A_3IM)*gdata.idx[0]);

	return Bphi;

#elif (GEOM == SPHERICAL) 

	REAL r_sph = gdata.get_x(ii    );
	REAL theta = gdata.get_y(jj+0.5);

	REAL rP = gdata.get_x(ii+0.5);
	REAL rM = gdata.get_x(ii-0.5);

	REAL Btheta = ((A_1KP    - A_1KM   )*gdata.idx[2]/sin(theta) -
	               (A_3IP*rP - A_3IM*rM)*gdata.idx[0])/r_sph;

	return Btheta;
	
#endif
#endif
}


REAL AlfvenWingsSpherical::B2(Data &gdata, int ii, int jj, int kk) {

#if (GEOM == CARTESIAN)

	REAL pos[3] = {ii+0.5,1.*jj,kk+0.5};
	REAL Bz = (Ay(gdata, pos) - Ay(gdata, pos))*gdata.idx[0];

	pos[0] = 1.*ii; pos[1] = jj+0.5; pos[2] = kk+0.5;
	Bz -= (Ax(gdata, pos) - Ax(gdata, pos))*gdata.idx[1];
	return Bz;

#else 

	REAL pos[3] = {ii+0.5,1.*jj,kk+0.5};
	REAL A_2IP = Get_x2Vec(gdata, Ax(gdata, pos), Ay(gdata, pos),
	                       Az(gdata, pos), pos);

	pos[0] = ii-0.5; pos[1] = 1.*jj; pos[2] = kk+0.5;
	REAL A_2IM = Get_x2Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);

	pos[0] = 1.*ii; pos[1] = jj+0.5; pos[2] = kk+0.5;
	REAL A_1JP = Get_x1Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);
	pos[0] = 1.*ii; pos[1] = jj-0.5; pos[2] = kk+0.5;
	REAL A_1JM = Get_x1Vec(gdata, Ax(gdata, pos), Ay(gdata, pos), 
	                       Az(gdata, pos), pos);

	REAL r_cyl = gdata.get_x(ii);

	REAL rP = gdata.get_x(ii+0.5);
	REAL rM = gdata.get_x(ii-0.5);

	REAL B3 = ((A_2IP*rP - A_2IM*rM)*gdata.idx[0] -
	           (A_1JP    - A_1JM   )*gdata.idx[2])/r_cyl;

#if (GEOM == CYLINDRICAL)

	REAL Bz = B3;
	return Bz;

#elif (GEOM == SPHERICAL) 

	REAL Bphi = B3;
	return Bphi;
	
#endif
#endif	
}



void AlfvenWingsSpherical::bc_User(Data &gdata, NumMatrix<REAL,3> &omb,
			int dir, int top, int q, int rim) 
{
  /*
    Include correct extrapolation for density and velocity for the
    initial conditions:
   */


  if(dir == 0) {
    int kmin = -1;
    if(q == 4) {
      kmin = -2;
    }
    if(top == 0) {
	    for (int k = -rim; k <= gdata.mx[2]+rim; ++k){
		    for (int j = -rim; j <= gdata.mx[1]+rim; j++){
			    for (int i = kmin; i >= -rim; --i){
				    if(q >= 0 && q < 8) {
					    omb(i,j,k) = GetConstFields(gdata, q,i,j,k);
				    } else {
					    omb(i,j,k) = omb(i+1,j,k);
				    }
			    }
		    }
      }
    } else {
      for (int k = -rim; k <= gdata.mx[2]+rim; ++k){
	      for (int j = -rim; j <=gdata. mx[1]+rim; ++j){
		      for (int i = 1; i <= rim; ++i){
			      if(q >= 0 && q < 8) {
				      omb(gdata.mx[0]+i,j,k) = GetConstFields(gdata, q,
				                                              gdata.mx[0]+i,j,k);
			      } else {
				      omb(i,j,k) = omb(gdata.mx[0]+i-1,j,k);
			      }	  
		      }
	      }
      }
    }

  } else if (dir == 1) {

	  int jmin = -1;
	  if(q == 5) {
      jmin = -2;
    }


    if(top == 0) {
      for (int k = -rim; k <= gdata.mx[2]+rim; k++){
	      for (int i = -rim; i <= gdata.mx[0]+rim; i++){
		      for (int j = jmin; j >= -rim; --j){
			      if(q >= 0 && q < 8) {
				      omb(i,j,k) = GetConstFields(gdata, q,i,j,k);
			      } else {
				      omb(i,j,k) = omb(i,j+1,k);
			      }
		      }
	      }
      }
    } else {
      for (int k = -rim; k <= gdata.mx[2]+rim; k++){
	      for (int i = -rim; i <= gdata.mx[0]+rim; i++){
		      for (int j = 1; j <= rim; j++){
			      if(q >= 0 && q < 8) {
				      omb(i,gdata.mx[1]+j,k) = GetConstFields(gdata, q,i,
				                                              gdata.mx[1]+j,k);
			      } else {
				      omb(i,gdata.mx[1]+j,k) = omb(i,gdata.mx[1]+j-1,k);
			      }
		      }
	      }
      }
    }

  } else if (dir == 2) {

	  int kmin = -1;
    if(q == 6) {
      kmin = -2;
    }

    if(top == 0) {
      for (int j = -rim; j <= gdata.mx[1]+rim; j++){
        for (int i = -rim; i <= gdata.mx[0]+rim; i++){
          for (int k = kmin; k >= -rim; --k){
            if(q >= 0 && q < 8) {
              omb(i,j,k) = GetConstFields(gdata, q,i,j,k);
            } else {
              omb(i,j,k) = omb(i,j,k+1);
            }
          }
        }
      }
    } else {
	    for (int j = -rim; j <= gdata.mx[1]+rim; ++j){
        for (int i = -rim; i <= gdata.mx[0]+rim; ++i){
          for (int k = 1; k <= rim; ++k){
            if(q >= 0 && q < 8) {
              omb(i,j,gdata.mx[2]+k) = GetConstFields(gdata, q,i,j,
                                                      gdata.mx[2]+k);
            } else {
              omb(i,j,gdata.mx[2]+k) = omb(i,j,gdata.mx[2]+k-1);
            }
          }
        }
	    }
    }

  } else {
    cerr << " Only applicable for x-direction " << endl;
    exit(414);
  }

}



#endif
