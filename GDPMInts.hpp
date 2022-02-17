#ifndef UNOMOL_GDPMINTS_HPP_
#define UNOMOL_GDPMINTS_HPP_
#include <iostream>
#include <cmath>
#include "Basis.hpp"
#include "AuxFunctions.hpp"
#include "Rys.hpp"
#include "Util.hpp"
#include "Structs.hpp"
using namespace std;

namespace unomol {

void calc_gdpm_ints(
    const ShellPairData& sq,
    double *svals,
    const AuxFunctions& aux,Rys& rys,const double *p);

void GDPMInts(const Basis& bas,double* Hmat);

}
#endif
