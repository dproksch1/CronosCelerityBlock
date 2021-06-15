#include "ConstrainedTransport.H"

void ConstrainedTransport::bc_EMF(Data &gdata, ProblemType &Problem,
                                  fields_2D &fieldsEmf, int dir) {
	// So far this routine is only used to apply some kind of boundary
	// conditions to the electric field, when the polar axis is
	// included in the simulations

	// Treatment of axis for cylindrical coordinates:
#if (GEOM == CYLINDRICAL)
	if(dir == 2) {
		if(gdata.get_singularity_treatment(0) == 1) { // at lower r-bound
			REAL emf_ave(0.);
			int nums(0);
#ifdef parallel
			if(gdata.coords[0] == 0) {
#endif
				// For every value of z sum over all values of phi at r=0
				for (int j = 0; j <= gdata.mx[1]; ++j){
					emf_ave += fieldsEmf.emf(0,j);
					nums++;
				}
#ifdef parallel
			}
			gdata.MpiSum(2, emf_ave); // Some only for ranks with same z
			gdata.MpiSum(2, nums);    // Some only for ranks with same z
			
			// gdata.MpiSum(emf_ave);
			// gdata.MpiSum(nums);
#endif
			emf_ave /= (1.*nums);
			// Apply average Ez instead of computed one
			for (int j = 0; j <= gdata.mx[1]+1; ++j){
				fieldsEmf.emf(0,j) = emf_ave;
			}
		} else if (gdata.get_singularity_treatment(0) == 2) { // at lower r 
#ifdef parallel
			// Only for ranks that are not located at the origin:
			REAL emf_ave(0.);
			int nums(0);
			gdata.MpiSum(2, emf_ave); // Some only for ranks with same z
			gdata.MpiSum(2, nums);    // Some only for ranks with same z
			
#endif
		}
	}
#endif


	// Treatment for spherical coordinates
#if (GEOM == SPHERICAL)

	// do axis boundary first
	if(dir == 0) {
		if(gdata.get_singularity_treatment(2) > 0) { // at lower theta-bound

			REAL emf_ave(0.);
			int nums(0);
#ifdef parallel
			if(gdata.coords[1] == 0) { // Check if rank at lower bound
#endif
				// Compute average over phi-direction
				for (int iphi = 0; iphi <= gdata.mx[2]; ++iphi){
					emf_ave += fieldsEmf.emf(0,iphi);
					nums++;
				}
#ifdef parallel
			}
			gdata.MpiSum(0, emf_ave); // Sum only for ranks with same r
			gdata.MpiSum(0, nums);    // Sum only for ranks with same r

			if(gdata.coords[1] == 0) { // Check if rank at lower bound
#endif	
				emf_ave /= (1.*nums);

				// Assign average Er instead of computed one
				for (int iphi = 0; iphi <= gdata.mx[2]+1; ++iphi){
					fieldsEmf.emf(0,iphi) = emf_ave;
				}
#ifdef parallel
			}
#endif
		
		}



		if(gdata.get_singularity_treatment(3) > 0) { // at upper theta-bound

			REAL emf_ave(0.);
			int nums(0);
#ifdef parallel
			// Check if rank at upper bound
			if(gdata.coords[1] == gdata.nproc[1]-1) {
#endif
				// Compute average over phi-direction
				for (int iphi = 0; iphi <= gdata.mx[2]; ++iphi){
					emf_ave += fieldsEmf.emf(gdata.mx[1]+1,iphi);
					nums++;
				}
#ifdef parallel
			}
			gdata.MpiSum(0, emf_ave); // Sum only for ranks with same r
			gdata.MpiSum(0, nums);    // Sum only for ranks with same r
			
			// Check if rank at upper bound
			if(gdata.coords[1] == gdata.nproc[1]-1) {
#endif	
				emf_ave /= (1.*nums);

				// Assign average Er instead of computed one
				for (int iphi = 0; iphi <= gdata.mx[2]+1; ++iphi){
					fieldsEmf.emf(gdata.mx[1]+1,iphi) = emf_ave;
				}
#ifdef parallel
			}
#endif
		}
	}


	// That one seems just wrong:
	if(false) {
	// do lower boundary first
	if(dir == 0) {
		if(gdata.get_singularity_treatment(0) == 1) { // at lower r-bound
			REAL emf_ave(0.);
			int nums(0);
#ifdef parallel
			if(gdata.coords[0] == 0) {
#endif
				// For every value of r sum over all values of phi at theta=0
				for (int k = 0; k <= gdata.mx[2]; ++k){
					emf_ave += fieldsEmf.emf(0,k);
					nums++;
				}
#ifdef parallel
			}
			gdata.MpiSum(0, emf_ave); // Some only for ranks with same r
			gdata.MpiSum(0, nums);    // Some only for ranks with same r
			
#endif
			emf_ave /= (1.*nums);
			// Apply average Ez instead of computed one
			for (int j = 0; j <= gdata.mx[1]+1; ++j){
				fieldsEmf.emf(0,j) = emf_ave;
			}
		} else if (gdata.get_singularity_treatment(0) == 2) { // at lower r 
#ifdef parallel
			// Only for ranks that are not located at the origin:
			REAL emf_ave(0.);
			int nums(0);
			gdata.MpiSum(0, emf_ave); // Some only for ranks with same r
			gdata.MpiSum(0, nums);    // Some only for ranks with same r
			
#endif
		}
	}
	}
#endif
	



#ifdef parallel
	MPI_Barrier(gdata.comm3d);
#endif
}



void ConstrainedTransport::bc_noms(Data &gdata, ProblemType &Problem, int dir) {
	// Set nom to zero at coordinate axes, when integrating the magnetic field
	// Currently, corrections is applied twice - because the error is produced in two directions
	if(!IntegrateA) {
#if (GEOM == CYLINDRICAL)
		if(dir==1 || dir==2) {
			if(gdata.get_singularity_treatment(0) == 1) { // at lower r-bound
				for (int j = -1; j <= gdata.mx[1]+1; ++j){
					for (int k = -1; k <= gdata.mx[2]; ++k){
						gdata.nom[q_Bx](-1,j,k) = 0.;
					}
				}
			}
		}
#elif(GEOM == IntegrateA)
		if(dir==0 || dir==2) {
			// fix axis for theta=0
			if(gdata.get_singularity_treatment(2) > 0) {
				for (int i_r = 0; i_r <= gdata.mx[0]; ++i_r){
					for (int i_phi = 0; i_phi <= gdata.mx[2]; ++i_phi){
						gdata.nom[q_By](i_r,-1,i_phi) = 0.;
					}
				}
			}
			// fix axis for theta=pi
			if(gdata.get_singularity_treatment(3) > 0) {
				for (int i_r = 0; i_r <= gdata.mx[0]; ++i_r){
					for (int i_phi = 0; i_phi <= gdata.mx[2]; ++i_phi){
						gdata.nom[q_By](i_r,gdata.mx[1],i_phi) = 0.;
					}
				}
			}
		}
#endif
	}
}
