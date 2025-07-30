#pragma once
#include <iostream>
#include <cmath>
#include "Basis.hpp"
#include "AuxFunctions.hpp"
#include "MD_Dfunction.hpp"
#include "Util.hpp"
#include "Structs.hpp"
using namespace std;

namespace unomol {
void AngMomInts(const Basis& bas);
void calculateAngMomentum(const double * pmat,size_t no2);
void calculateAngMomentum(const double * pmatA, const double * pmatB, size_t no2);
}
