#ifndef CRONOS_CONSTRAINED_TRANSPORT_H
#define CRONOS_CONSTRAINED_TRANSPORT_H = 1

#include "vector.H"
#include "fields_local.H"
#include "DissipationMHD.H"
#include "reconst.H"
#include "reconst_2nd.H"
#include "save_data.H"
#include "transformations.H"

class ConstrainedTransport {
public:
	ConstrainedTransport(const Data &, int dir, bool IntegrateA, int i_magFluid=0);
	virtual ~ConstrainedTransport();
//#if (USE_COROTATION == CRONOS_ON)
//	void get_NumEMFCT(Data &, ProblemType &, Saves &, Transformations &Trafo,
//	                  NumMatrix<REAL, 3> nom [N_OMINT]);
//#else
	void get_NumEMFCT(Data &, ProblemType &, Saves &,
	                  NumMatrix<REAL, 3> nom[N_OMINT]);
//#endif
	void gen_EmfData(Data &);
protected:
	void get_ChangesEmf(Data &,
	                    NumMatrix<REAL, 3> nom [N_OMINT],
	                    int layer);
#if (USE_COROTATION == CRONOS_ON)
	void TransInertToCorot(Data &gdata, Transformations &Trafo, int layer);
#endif
	virtual void get_LocalArrays(const Data &, int layer);
	virtual void Reconstruct(const Data &){}
	virtual void get_PhysEmfs(const Data &, const ProblemType &,
	                          phys_fields_2D &, cronos::vector<REAL> &ipos){}
	virtual void get_vChar2D(const Data &, const ProblemType &,
	                         cronos::vector<REAL> &ipos){}
	virtual void get_NumEmf2D(const Data &, const phys_fields_2D &physValEmfLL,
	                          const phys_fields_2D &physValEmfLR,
	                          const phys_fields_2D &physValEmfRL,
	                          const phys_fields_2D &physValEmfRR,
	                          fields_2D &fieldsEmf) = 0;
	virtual void bc_EMF(Data &gdata, ProblemType &Problem,
	                    fields_2D &fieldsEmf, int dir);
	/** Apply boundary conditions to the noms after application of changes by
	 * constrained transport -> currently used to avoid NANs on z-axis in
	 * cylindrical and spherical coordinates */
	virtual void bc_noms(Data &gdata, ProblemType &Problem, int dir);

	phys_fields_2D *physValEmfLL, *physValEmfLR, *physValEmfRL, *physValEmfRR;
	Reconstruction2D *ReconstEmf;
	fields_2D *fieldsEmf;
	DissipationMHD *DissMHD;

	cronos::vector<REAL> ipos;

	int dir, dir0, dir1;
	bool IntegrateA;
	int q_sx, q_sy, q_sz;
	int q_Bx, q_By, q_Bz;
	//	bool with_axis;
};

#endif
