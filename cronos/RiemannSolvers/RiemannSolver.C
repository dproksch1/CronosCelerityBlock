#include "RiemannSolver.H"
#include <iostream>

using namespace std;

void RiemannSolver::compute_carbuncleFlag(Data &gdata) {
	//! Compute flag to avoid carbunce problem
	/*!
	 * To circumvent the carbuncle problem present in hllc and hlld
	 * we use the procedure suggested in Quirk (1992)
	 * Quirk suggests to use a flag to indicate the presence of a strong
	 * shock. Only in those cells, where such a shock is present, the
	 * Riemann solver is reduced to hll
	 * Parameter alpha, an estimate of the pressure gradient at a strong shock,
	 * is to be supplied by the user as 'alpha_carbuncle'
	 * RK 2014-12-16
	 * */

//	if(gdata.use_carbuncleFlag) {
//		double pLoc, pOth1, pOth2, ratio;
//		// Loop over whole grid:
//		for(int zk = -2; zk<=gdata.mx[2]+2; ++zk){
//			for(int jy = -2; jy<=gdata.mx[1]+2; ++jy){
//				for(int ix = -2; ix<=gdata.mx[0]+2; ++ix){
//					int local_flag = 0;
//					pLoc = gdata.pTherm(ix, jy, zk);
//
//					for(int dir=0; dir<3; ++dir) {
//						if(dir==0) {// x-direction
//							pOth1 = gdata.pTherm(ix+1, jy, zk);
//							pOth2 = gdata.pTherm(ix-1, jy, zk);
//						} else if (dir==1) {
//							pOth1 = gdata.pTherm(ix, jy+1, zk);
//							pOth2 = gdata.pTherm(ix, jy-1, zk);
//						} else {
//							pOth1 = gdata.pTherm(ix, jy, zk+1);
//							pOth2 = gdata.pTherm(ix, jy, zk-1);
//						}
//						ratio = abs(pLoc-pOth1)/std::min(pLoc, pOth1);
//						if(ratio > alpha_carbuncle) {
//							local_flag = 1;
//						}
//						ratio = abs(pLoc-pOth2)/std::min(pLoc, pOth2);
//						if(ratio > alpha_carbuncle) {
//							local_flag = 1;
//						}
//					}
//
//					gdata.carbuncleFlag(ix, jy, zk) = local_flag;
//
//				}
//			}
//		}
//		cout << " The flag ";
//		cout << gdata.carbuncleFlag(23,182,0) << " ";
//		cout << gdata.carbuncleFlag(23,183,0) << " ";
//		cout << gdata.carbuncleFlag(23,184,0) << " ";
//		cout << gdata.carbuncleFlag(23,185,0) << " ";
//		cout << endl;
//	}
	if(gdata.use_carbuncleFlag) {

		gdata.carbuncleFlag.clear();

		double pLoc, pOth, ratio;
		// Loop over whole grid:
		for(int zk = -2; zk<=gdata.mx[2]+2; ++zk){
			for(int jy = -2; jy<=gdata.mx[1]+2; ++jy){
				for(int ix = -2; ix<=gdata.mx[0]+2; ++ix){
					int local_flag = 0;
					pLoc = gdata.pTherm(ix, jy, zk);

					for(int dir=0; dir<3; ++dir) {
						if(dir==0) {// x-direction
							pOth = gdata.pTherm(ix+1, jy, zk);
						} else if (dir==1) {
							pOth = gdata.pTherm(ix, jy+1, zk);
						} else {
							pOth = gdata.pTherm(ix, jy, zk+1);
						}
						ratio = abs(pLoc-pOth)/std::min(pLoc, pOth);
//						if(ratio>0.1) {
//							cout << " Da " << ratio << " " << ix << " " << jy << endl;
//						}
						if(ratio > alpha_carbuncle) {
							local_flag = 1;
//							if(zk==0) {
//								cout << " Flag at " << ix << " " << jy << " " << zk << " " << dir <<  endl;
//							}
							if(dir==0) {
								gdata.carbuncleFlag(ix, jy, zk) = 1;
								gdata.carbuncleFlag(ix+1, jy, zk) = 1;
							} else if(dir==1) {
								gdata.carbuncleFlag(ix, jy, zk) = 1;
								gdata.carbuncleFlag(ix, jy, zk+1) = 1;
							} else {
								gdata.carbuncleFlag(ix, jy, zk) = 1;
								gdata.carbuncleFlag(ix, jy, zk+1) = 1;
							}
						}
					}

//					gdata.carbuncleFlag(ix, jy, zk) = local_flag;

				}
			}
		}
//		cout << " The flag ";
//		cout << gdata.carbuncleFlag(23,182,0) << " ";
//		cout << gdata.carbuncleFlag(23,183,0) << " ";
//		cout << gdata.carbuncleFlag(23,184,0) << " ";
//		cout << gdata.carbuncleFlag(23,185,0) << " ";
//	cout << endl;
	}

}


void RiemannSolver::set_verbosity(int _verb) {
	verbosity = _verb;
}
