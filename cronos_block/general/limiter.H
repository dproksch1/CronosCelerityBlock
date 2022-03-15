#ifndef CRONOS_LIMITER_H
#define CRONOS_LIMITER_H = 1

#include <string>
#include <variant>

#include "constants.H"
#include "minmod.H"

class limiter {
public:
	limiter(int substep = 0) {
		int limiterID = -1;
		//if (value_exists("Limiter_" + ToString(substep)))
		//{
		//	limiterID = static_cast<int>(value("Limiter_" + ToString(substep)));
		//}
		//else
		//{
		//	limiterID = static_cast<int>(value("Limiter"));
		//}
		createLimiter(limiterID);
	}

	limiter(const limiter&) = default;

	double compute(double& up, double& uc, double& um) const {
		return Limiter.compute(up, uc, um);
		// ugly workaround because there is no noexcept version of std::variant...
		//if (auto p = std::get_if<minmod>(&Limiter)) {
		//	return p->compute(up, uc, um);
		//} else if (auto p = std::get_if<vanLeer>(&Limiter)) {
		//	return p->compute(up, uc, um);
		//} else if (auto p = std::get_if<Superbee>(&Limiter)) {
		//	return p->compute(up, uc, um);
		//}
		//return 0.0;
		//return std::visit([&](auto&& arg) { return arg.compute(up, uc, um); }, Limiter);
	}

	std::string get_Name() {
		return name;
	}

protected:
	//void set_limiter(int);

	void createLimiter(int ch_limit) {
		Limiter = minmod();
		//if (ch_limit == 1) {
		//	Limiter = minmod();
		//	name = "minmod";
		//}
		//else if (ch_limit == 2) {
		//	Limiter = vanLeer();
		//	name = "van Leer";
		//}
		//else if (ch_limit == 3) {
		//	Limiter = Superbee();
		//	name = "Superbee";
		//}
		//else {
		//	//cout << " No such limiter available " << endl;
		//	//cout << " Choices: " << endl;
		//	//cout << "    minmod -->  (1) " << endl;
		//	//cout << "  van Leer -->  (2) " << endl;
		//	//cout << "  Superbee -->  (3) " << endl;
		//}
	}

	//std::variant<minmod, vanLeer, Superbee> Limiter;
	minmod Limiter;

	std::string name;
};


#endif
