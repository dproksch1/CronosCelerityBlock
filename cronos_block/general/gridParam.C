#include "gridParam.H"
#include <iostream>

GridFunction_Exp::GridFunction_Exp(double domain_begin, double domain_end,
                                   double domain_len) {
	this->domain_begin = domain_begin;
	this->domain_end = domain_end;
	this->domain_len = domain_len;
}

double GridFunction_Exp::get_gridFunc(double ratio) {
	return domain_begin/domain_len*
		(pow(domain_end/domain_begin,ratio) - 1.);
}


GridFunction_Log::GridFunction_Log(double domain_begin, double domain_end) {
	this->domain_begin = domain_begin;
	this->domain_end = domain_end;
}

double GridFunction_Log::get_gridFunc(double ratio) {
	std::cout << domain_begin*pow(domain_end/domain_begin,ratio) << std::endl;
	return domain_begin*pow(domain_end/domain_begin,ratio);
}


GridFunction_Ratio::GridFunction_Ratio(double ref, double min_ratio, double max_ratio) {
	_ref = ref;
	_min_ratio = min_ratio;
	_max_ratio = max_ratio;
}

double GridFunction_Ratio::get_gridFunc(double ratio) {
	return _min_ratio*pow(_max_ratio/_min_ratio,ratio)*_ref + _ref;
}


GridFunction_Percent::GridFunction_Percent(double domain_beg, double domain_end, double increase, int num) {

	_domain_beg = domain_beg;
	_increase = increase;
	_num = num;

	_delta = (domain_end - domain_beg)*(1-increase)/(1.-pow(increase,num));

}

double GridFunction_Percent::get_gridFunc(double ratio) {
	double iVal = ratio*_num;
//	if(iVal <20 || iVal > 320) {
//		std::cout << "grid: " << iVal << " " << ratio << " " << _num << " ";
//		std::cout << _domain_beg + _delta*(1-pow(_increase, iVal))/(1.-_increase) << " ";
//		std::cout << _delta<< " ";
//		std::cout << _delta/_domain_beg << " ";
//		std::cout << std::endl;
//	}
	return _domain_beg + _delta*(1-pow(_increase, iVal))/(1.-_increase);
}

GridFunction_PercentVar::GridFunction_PercentVar(double domain_beg, double delta, double increase, int num) {
	/*!
	 * Constructor of non-linear grid similar to the one described in ud-Doula & Owocki (2002) on page 4
	 * @param[in] domain_beg  beginning of numerical domain - left face of cell 0
	 * @param[in] delta       size of cell 0
	 * @param[in] increase    change of cell size from cell i to cell i+1
	 * @param[in] num         total number of grid points
	 */
	_domain_beg = domain_beg;
	_increase = increase;
	_num = num;
	_delta = delta;

}

double GridFunction_PercentVar::get_gridFunc(double ratio) {
	double iVal = ratio*_num;
//	if(iVal <20) {
//		std::cout << "grid: " << iVal << " " << ratio << " " << _num << " ";
//		std::cout << _domain_beg + _delta*(1-pow(_increase, iVal))/(1.-_increase) << " ";
//		std::cout << _delta<< " ";
//		std::cout << _delta/_domain_beg << " ";
//		std::cout << std::endl;
//	}
	return _domain_beg + _delta*(1-pow(_increase, iVal))/(1.-_increase);
}
