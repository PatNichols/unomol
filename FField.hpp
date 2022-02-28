#ifndef UNOMOL_FFIELD_hpp_
#define UNOMOL_FFIELD_hpp_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "Basis.hpp"
#include "AuxFunctions.hpp"
#include "MD_Dfunction.hpp"
#include "MD_Rfunction.hpp"
#include "Util.hpp"
#include "Structs.hpp"
#include "putils_c.h"
using namespace std;

namespace unomol {
void FiniteFieldMatrix( const Basis& basis, double *Ham,
                        const double *Efield);

}

#endif


