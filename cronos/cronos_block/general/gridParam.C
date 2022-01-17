#include "gridParam.H"
#include <iostream>

GridFunction_Exp::GridFunction_Exp(REAL domain_begin, REAL domain_end,
                                   REAL domain_len) {
	this->domain_begin = domain_begin;
	this->domain_end = domain_end;
	this->domain_len = domain_len;
}

REAL GridFunction_Exp::get_gridFunc(REAL ratio) {
	return domain_begin/domain_len*
		(pow(domain_end/domain_begin,ratio) - 1.);
}


GridFunction_Log::GridFunction_Log(REAL domain_begin, REAL domain_end) {
	this->domain_begin = domain_begin;
	this->domain_end = domain_end;
}

REAL GridFunction_Log::get_gridFunc(REAL ratio) {
	std::cout << domain_begin*pow(domain_end/domain_begin,ratio) << std::endl;
	return domain_begin*pow(domain_end/domain_begin,ratio);
}


GridFunction_Ratio::GridFunction_Ratio(REAL ref, REAL min_ratio, REAL max_ratio) {
	_ref = ref;
	_min_ratio = min_ratio;
	_max_ratio = max_ratio;
}

REAL GridFunction_Ratio::get_gridFunc(REAL ratio) {
	return _min_ratio*pow(_max_ratio/_min_ratio,ratio)*_ref + _ref;
}


GridFunction_Percent::GridFunction_Percent(REAL domain_beg, REAL domain_end, REAL increase, int num) {

	_domain_beg = domain_beg;
	_increase = increase;
	_num = num;

	_delta = (domain_end - domain_beg)*(1-increase)/(1.-pow(increase,num));

}

REAL GridFunction_Percent::get_gridFunc(REAL ratio) {
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

GridFunction_PercentVar::GridFunction_PercentVar(REAL domain_beg, REAL delta, REAL increase, int num) {
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

REAL GridFunction_PercentVar::get_gridFunc(REAL ratio) {
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
