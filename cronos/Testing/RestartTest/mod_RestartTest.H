#ifndef CRONOS_RESTART_TEST
#define CRONOS_RESTART_TEST

#include "problem.H"
#include "randgen.H"

// using namespace std;

class RestartTest : public ProblemType {
 public:
  RestartTest(Data &);
  virtual void init_fields(Data &, int[], int[]);

 private:
  RandomNumberGenerator RandGen;
};

RestartTest::RestartTest(Data &gdata) : ProblemType(gdata) {}

void RestartTest::init_fields(Data &gdata, int ibeg[3], int iend[3]) {
  for (int iz = ibeg[2]; iz <= iend[2]; ++iz) {
    for (int iy = ibeg[1]; iy <= iend[1]; ++iy) {
      for (int ix = ibeg[0]; ix <= iend[0]; ++ix) {
        gdata.om[q_rho](ix, iy, iz) = RandGen.getRand();
        gdata.om[q_sx](ix, iy, iz) = RandGen.getRand();
        gdata.om[q_sy](ix, iy, iz) = RandGen.getRand();
        gdata.om[q_sz](ix, iy, iz) = RandGen.getRand();
        gdata.om[q_Eges](ix, iy, iz) = RandGen.getRand();
      }
    }
  }

#if (FLUID_TYPE == CRONOS_MHD)
  for (int iz = ibeg[2]; iz <= iend[2]; ++iz) {
    for (int iy = ibeg[1]; iy <= iend[1]; ++iy) {
      for (int ix = ibeg[0]; ix <= iend[0]; ++ix) {
        gdata.om[q_Bx](ix, iy, iz) = RandGen.getRand();
        gdata.om[q_By](ix, iy, iz) = RandGen.getRand();
        gdata.om[q_Bz](ix, iy, iz) = RandGen.getRand();
      }
    }
  }

  gdata.om[q_Bx].rename("A_x");
  gdata.om[q_By].rename("A_y");
  gdata.om[q_Bz].rename("A_z");
#endif
}

#endif
