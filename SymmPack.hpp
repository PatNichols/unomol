#ifndef UNOMOL_SymmPack_HPP
#define UNOMOL_SymmPack_HPP
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "Util.hpp"
using namespace std;

namespace unomol {
namespace SymmPack {

double TraceSymmPackProduct(
    double a[],double b[],int n) noexcept;

double SymmPackDiffNorm(
    double a[],double b[],int n) noexcept;

void rsp(const int n,
         double* __restrict__ a,
         double* __restrict__ z,
         double* __restrict__ evals,
         double* __restrict__ tmp) noexcept;


void sp_trans (int n, double *a, double *z, double *tmp) noexcept;

}  //end namespace
} //end namespace
#endif


