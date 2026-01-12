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

void GDPMInts(const Basis& bas, double* Hmat);

}
#endif
