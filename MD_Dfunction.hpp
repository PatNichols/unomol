#ifndef UNOMOL_MD_DFUNCTION_HPP
#define UNOMOL_MD_DFUNCTION_HPP
#include <iostream>
#include <cmath>
#include "Util.hpp"
using namespace std;

namespace unomol {

class MD_Dfunction {

  private:
    int asize1;
    double*** dtx;

  public:
    MD_Dfunction(int lmax) {
        asize1=lmax+1;
        int asize2=2*lmax+1;
        dtx=new_tensor3< double >(asize1,asize1,asize2);
    }

    ~MD_Dfunction() {
        delete_tensor3<double>(dtx,asize1,asize1);
    }

    constexpr const double* getRow(int i,int j) const noexcept {
        return dtx[i][j];
    }

    constexpr double getValue(int i,int j,int k) const noexcept {
        return dtx[i][j][k];
    }

    constexpr const double * operator()(int i,int j) const noexcept { return dtx[i][j];}

    constexpr double operator()(int i,int j,int k) const noexcept { return dtx[i][j][k];}

    constexpr void eval(double abi,double ax,double bx,int l1,int l2) noexcept {
        int ltot=l1+l2;
        if (ltot==0) {
            dtx[0][0][0]=1.0;
            return;
        }
        for (int i=0; i<=l1; i++) {
            for (int j=0; j<=l2; j++) {
                for (int k=0; k<=ltot; k++) dtx[i][j][k]=0.0;
            }
        }
        dtx[0][0][0]=1.0;
        for (int i=1; i<=l2; i++) {
            int im1=i-1;
            dtx[0][i][0]=bx*dtx[0][im1][0]+dtx[0][im1][1];
            for (int n=1; n<i; n++) {
                dtx[0][i][n]=abi*dtx[0][im1][n-1]+bx*dtx[0][im1][n]+
                             (n+1)*dtx[0][im1][n+1];
            }
            dtx[0][i][i]=abi*dtx[0][im1][im1];
        }
        for (int j=1; j<=l1; j++) {
            for (int i=0; i<=l2; i++) {
                int jm1=j-1;
                int ipj=i+j;
                dtx[j][i][0]=ax*dtx[jm1][i][0]+dtx[jm1][i][1];
                for (int n=1; n<ipj; n++) {
                    dtx[j][i][n]=abi*dtx[jm1][i][n-1]+ax*dtx[jm1][i][n]+
                                 (n+1)*dtx[jm1][i][n+1];
                }
                dtx[j][i][ipj]=abi*dtx[jm1][i][ipj-1];
            }
        }
    }
};

}
#endif
