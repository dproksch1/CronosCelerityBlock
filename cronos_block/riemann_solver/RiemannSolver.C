#include "RiemannSolver.H"
#include <iostream>

using namespace std;

void RiemannSolver::compute_carbuncleFlag(Data &gdata) {
	
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
				}
			}
		}
	}
}


void RiemannSolver::set_verbosity(int _verb) {
	verbosity = _verb;
}
