#ifndef UNOMOL_MOMENTS_HPP
#define UNOMOL_MOMENTS_HPP
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

void calc_moments (MomInts * mvals,
                   ShellPairData & sp,
                   const AuxFunctions & aux,
                   MD_Dfunction & dx, MD_Dfunction & dy, MD_Dfunction & dz);

void MomentInts (const Basis & basis );

void AnalyzeMoments (const double *Pmat, const Center * center, int ncen, int no2);

void AnalyzeMoments (const double *PmatA, const double* PmatB,
                     const Center * center, int ncen, int no2);

void AnalyzeMOMoments (double *Cmat,int no2, int norb);

void AnalyzeMOMoments (double *CmatA,double *CmatB,int no2, int norb);

}
#endif
