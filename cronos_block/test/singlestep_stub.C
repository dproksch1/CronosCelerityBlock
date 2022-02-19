#include "gridgen.H"
#include <iomanip>
#include <cmath>
#include <sys/time.h>
#include <vector>
//#include "DissipationMHD.H"
#include "reconst_block.H"
#include <CL/sycl.hpp>

using namespace std;

double HyperbolicSolver::singlestep(Data &gdata, gridFunc &gfunc,
                                  ProblemType &Problem, int n, Queue&)
{
    //do nothing
    return 0.0;
}