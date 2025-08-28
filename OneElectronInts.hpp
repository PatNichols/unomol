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

void  OneElectronInts(const Basis& bas,double* Smat,
                      double* Tmat, double* Hmat);

}
#endif
