#ifndef UNOMOL_AUX_FUNCTIONS_hpp
#define UNOMOL_AUX_FUNCTIONS_hpp
#include <cstdlib>
#include <cmath>
#include <cstring>

namespace unomol {

/////////////////////////////////////////////////////
// Auxillary Functions used for doing calculations //
//  most of these derive from the need to use      //
//  shells instead of true orbitals                //
class AuxFunctions {
  private:
    int maxlst;
    int maxl;
    int *nlst;
    int ***lxyz;
    double **nfact;
  public:
    AuxFunctions(int maxl=2) {
        int mlp1=maxl+1;
        nlst=new int[mlp1];
        for (int i=0; i<=maxl; i++) {
            nlst[i]=((i+1)*(i+2))/2;
        }
        maxlst=nlst[maxl];
        lxyz=new int**[mlp1];
        for (int i=0; i<=maxl; i++) {
            int nls=nlst[i];
            lxyz[i]=new int*[nls];
            for (int lx=0; lx<nls; ++lx) lxyz[i][lx]=new int[3];
            int k=0;
            for (int lx=i; lx>=0; lx--) {
                for (int ly=i-lx; ly>=0; ly--) {
                    lxyz[i][k][0]=lx;
                    lxyz[i][k][1]=ly;
                    lxyz[i][k][2]=i-lx-ly;
                    k++;
                }
            }
        }
        nfact=new double*[mlp1];
        for (int i=0; i<=maxl; i++) {
            int nls=nlst[i];
            nfact[i]=new double[nls];
        }
        double * dfact = new double[mlp1];
        dfact[0] = 1.;
        double dx = 1.;
        for (int i=1; i<=maxl; ++i) {
            dfact[i] = dfact[i-1] * dx;
            dx *= (2*i+1);
        }
        for (int i=0; i<=maxl; i++) {
            int nls=nlst[i];
            for (int j=0; j<nls; j++) {
                double dx=dfact[lxyz[i][j][0]];
                double dy=dfact[lxyz[i][j][1]];
                double dz=dfact[lxyz[i][j][2]];
                nfact[i][j]=1.0/sqrt(dx*dy*dz);
            }
        }
        delete [] dfact;
    }

    ~AuxFunctions() {
        int maxlp1 = maxl + 1;
        for (int l=maxlp1; l;) {
            --l;
            delete [] nfact[l];
        }
        delete [] nfact;
        for (int l=maxlp1; l;) {
            --l;
            int nls = nlst[l];
            for (int is=(nls+1); is;) {
                --is;
                delete [] lxyz[l][is];
            }
            delete [] lxyz[l];
        }
        delete [] lxyz;
        delete [] nlst;
    }

    constexpr int maxLstates() const noexcept {
        return maxlst;
    }
    // number of L vectors for lvalue lv //
    constexpr int number_of_lstates(int lv) const noexcept {
        return nlst[lv];
    }
    // the L vector for lvalue lv and l state ls //
    constexpr const int* l_vector(int lv,int ls) const noexcept {
        return lxyz[lv][ls];
    }
    // tensor for storing the L vectors for lvalue lv and l state ls //
    constexpr int Lxyz(int lv,int ls,int i) const noexcept {
        return lxyz[lv][ls][i];
    }
    // vector of L dependent normalization factors //
    constexpr const double* norm_factor_vec(int lv) const noexcept {
        return nfact[lv];
    }
    // L vector dependent normalization factor //
    constexpr double normalization_factor(int lv,int ls) const noexcept {
        return nfact[lv][ls];
    }
};

}
#endif
