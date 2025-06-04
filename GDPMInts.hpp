#ifndef UNOMOL_GDPMINTS_HPP_
#define UNOMOL_GDPMINTS_HPP_
#include <iostream>
#include <cmath>
#include "Basis.hpp"
#include "AuxFunctions.hpp"
#include "MD_Dfunction.hpp"
#include "MD_Rfunction.hpp"
#include "Util.hpp"
#include "Structs.hpp"

using namespace std;

namespace unomol {

double calc_dpm_nucrep(const Center * center, int ncen, int skipcen);

void GDPMInts(const Basis& bas,double* Hmat);

}
#endif
