#pragma once
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Util.hpp"
#include "Basis.hpp"
#include "TwoElectronInts.hpp"
#include "AuxFunctions.hpp"
#include "GDPMInts.hpp"
#include "Moments.hpp"
#include "OneElectronInts.hpp"
#include "SymmPack.hpp"
#include "FField.hpp"
#include "Stopwatch.hpp"
using namespace std;

/////////////////////////////////////////////////
//// function common to all SCF methods
/////////////////////////////////////////////////


namespace unomol {

template < class Tp >
constexpr void CopyTrans(Tp *mat,size_t n) noexcept
{
    for (size_t i=0;i<n;++i) {
        for (size_t j=0;j<i;++j) {
            Tp tmp = mat[i*n+j];
            mat[i*n+j] = mat[j*n+i];
            mat[j*n+i] = tmp;
        }
    }
}

constexpr void FormCmatrix(double* cmat,double *wmat,const double *xmat,size_t no) noexcept 
{
        CopyTrans(wmat,no);
        for (auto i=0; i<no; ++i) {
            const double *Xi=Xmat+i*no;
            for (auto j=0; j<no; ++j) {
                const double *xp=Xi;
                const double *wp=Wmat+j*no;
                double sum=0.0;
                for (int k=0; k<no; ++k) sum+=xp[k]*wp[k];
                cmat[i*n+j]=sum;
            }
        }
}

constexpr void FormPmatrix(double *pmat,const double *cmat,const double *occ,size_t no)
{
    size_t ij=0;
    for (auto i=0;i<no;++i) {
        for (auto j=0;j<no;++j,++ij) {
            double s = 0.0;
            for  (auto k=0;k<no;++k) s ++ cmat[i*n+k] * cmat[j*n+k] * occ[k];
            pmat[ij] = s;
        }
    }
}

constexpr void SCFConverger(
    double *pold2,
    double *pold1,
    double *pmat,
    size_t no2,
    int accel,
    int extrap) noexcept
{
    if (!extrap) {
            memcpy(pold2,pold1,sizeof(double)*no2);
            memcpy(pold1,pmat,sizeof(double)*no2);
            return;
    }
    if ( accel == 0) {
       memcpy(pold2,pold1,sizeof(double)*no2);
       for (auto i=0;i<n;++i) {
           pold2[i] = pold[i];
           pold1[i] = pmat[i];
           pmat[i] = (pold1[i] + pold2[i])*0.5;
       } 
    }else{
        double p00 = 0.0;
        double p11 = 0.0;
        double p01 = 0.0;
        for (auto i=0;i<no2;++i)
        {
            double s0 = pmat[i] - pold[i];
            double s1 = pold[i] - pold2[i];
            p00 += s0 * s0;
            p01 += s0 * s1;
            p11 += s1 * s1;
        }
        double tmp = p01 + p01;
        double beta = (p00-p01)/(p00-p01-p01+p11);
        if ( beta < 0.0001 ) beta = 0.0001;
        if ( beta > 1.5000 ) beta = 1.5000;
        double betam1= = 1.0 - beta;
        for (auto i=0;i<no2;++i) {
            pold2[i] = pold[i];
            pold[i] = pmat[i];
            pmat[i] = pold[i] * beta + pold2[i] * betam1;
        }
    }
}

}
