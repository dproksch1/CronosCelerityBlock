#ifndef CRONOS_CT_LONDRILLO_DELZANNA_H
#define CRONOS_CT_LONDRILLO_DELZANNA_H = 1

#include "ConstrainedTransport.H"

class CTLondrilloDelZanna: public ConstrainedTransport {
public:
	CTLondrilloDelZanna(const Data &, int dir, bool IntegrateA, int i_magFluid=0);
	virtual ~CTLondrilloDelZanna(){}
protected:
	virtual void Reconstruct(const Data &);
	virtual void get_PhysEmfs(const Data &, const ProblemType &,
	                          phys_fields_2D &, cronos::vector<REAL> &ipos);
	virtual void get_vChar2D(const Data &, const ProblemType &,
	                         cronos::vector<REAL> &ipos);
	virtual void get_NumEmf2D(const Data &, const phys_fields_2D &physValEmfLL,
	                          const phys_fields_2D &physValEmfLR,
	                          const phys_fields_2D &physValEmfRL,
	                          const phys_fields_2D &physValEmfRR,
	                          fields_2D &fieldsEmf);
private:
	REAL veps;
};


#endif
