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

void RiemannSolver::compute_carbuncleFlag(Queue &queue, Data &gdata) {
	
	if(gdata.use_carbuncleFlag) {

		double alpha_c = alpha_carbuncle;
		auto range = gdata.omSYCL[0].get_range();

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

			celerity::accessor ptherm_acc{gdata.pThermSYCL[0], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};	
			celerity::accessor carbuncleFlag_acc{gdata.carbuncleFlagSYCL[0], cgh, celerity::access::neighborhood{2,2,2}, celerity::write_only};

			cgh.parallel_for<class CarbuncleFlagKernel>(range, [=](celerity::item<3> item) {

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				double pLoc, pOth, ratio;

				int local_flag = 0;

				pLoc = ptherm_acc[ix][iy][iz];

				for(int dir=0; dir<3; ++dir) {
					if(dir==0) {// x-direction
						pOth = ptherm_acc[ix+1][iy][iz];
					} else if (dir==1) {
						pOth = ptherm_acc[ix][iy+1][iz];
					} else {
						pOth = ptherm_acc[ix][iy][iz+1];
					}
					ratio = abs(pLoc-pOth)/cl::sycl::min(pLoc, pOth);

					if(ratio > alpha_c) {
						if(dir==0) {
							carbuncleFlag_acc[ix][iy][iz] = 1;
							carbuncleFlag_acc[ix+1][iy][iz] = 1;
						} else if(dir==1) {
							carbuncleFlag_acc[ix][iy][iz] = 1;
							carbuncleFlag_acc[ix][iy][iz+1] = 1;
						} else {
							carbuncleFlag_acc[ix][iy][iz] = 1;
							carbuncleFlag_acc[ix][iy][iz+1] = 1;
						}
					}
				}
			});
		});
	}
}


void RiemannSolver::set_verbosity(int _verb) {
	verbosity = _verb;
}
