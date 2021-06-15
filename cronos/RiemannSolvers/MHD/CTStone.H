#ifndef CRONOS_CT_STONE_H
#define CRONOS_CT_STONE_H = 1

#include "ConstrainedTransport.H"

class CTStone: public ConstrainedTransport {
public:
	CTStone(const Data &, int dir, bool IntegrateA, int i_magFluid=0);
	virtual ~CTStone(){}
protected:
	virtual void get_NumEmf2D(const Data &, const phys_fields_2D &physValEmfLL,
	                          const phys_fields_2D &physValEmfLR,
	                          const phys_fields_2D &physValEmfRL,
	                          const phys_fields_2D &physValEmfRR,
	                          fields_2D &fieldsEmf);
	REAL eps_CentreCT;
};



#endif
