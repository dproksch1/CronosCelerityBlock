#include "PhysFluxes.H"
#include "pde_Transport.H"
#include "specific.H"




void HyperbolicSolver::set_UserPde(const Data &gdata) {

	/*
	  This routine can be modified to add different cases for example.
	 */

	PhysFluxUser = new PhysFluxUserTransport;


}
