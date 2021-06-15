#include "../Orbit.H"
#include "iostream"
#include <string>

using namespace std;

double computeOppositePhase(double supcPhase, double eccentricity) {
  const double eccAnomaly = ORBIT::PHASE::eccentricAnomaly(supcPhase, eccentricity);
  const double trueAnomaly = ORBIT::PHASE::trueAnomaly(supcPhase, eccentricity);
  const double trueAnomalyOpposite = trueAnomaly + M_PI;
  const double oppositePhase = ORBIT::PHASE::trueAnomaly2Phase(trueAnomalyOpposite, eccentricity);
  return oppositePhase;
}

template<typename T>
int EQ(string name, T a, T b, T tol = 0) {
  if (abs(a - b) <= tol * 0.5 * abs(a + b)) {
    cout << "OK: " << name << endl;
    return 0;
  } else {
    cout << "FAILED: " << name << endl;
    return 1;
  }
}

int main() {
  int ret = 0;

  ret += EQ("SUPC1", computeOppositePhase(0.08, 0.35), 0.76984, 1e-5);
  ret += EQ("SUPC2", computeOppositePhase(0.058, 0.35), 0.716153, 1e-5);
  ret += EQ("SUPC3", computeOppositePhase(0.06, 0.24), 0.649987, 1e-5);

  return ret;
}
