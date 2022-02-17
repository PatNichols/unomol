#ifndef UNOMOL_ONE_ELECTRON_INTS_HPP
#define UNOMOL_ONE_ELECTRON_INTS_HPP
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

void calc_one_electron_ints(
    const ShellPairData& sp,
    double svals[],double tvals[],double vvals[],
    const Center* center,int ncen,int skip,
    const AuxFunctions& aux,
    MD_Dfunction& dx,MD_Dfunction& dy,MD_Dfunction& dz,
    MD_Rfunction& r,double*** rsum);

void  OneElectronInts(const Basis& bas,double* Smat,
                      double* Tmat, double* Hmat);

}
#endif
