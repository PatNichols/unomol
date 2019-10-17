#ifndef _AUX_FUNCTIONS_CC_
#define _AUX_FUNCTIONS_CC_
#include <cstdlib>
#include <cmath>
#include <cstring>
using namespace std;
/////////////////////////////////////////////////////
// Auxillary Functions used for doing calculations //
//  most of these derive from the need to use      //
//  shells instead of true orbitals                //
class AuxFunctions
{
private:
    int maxlst;
    int *nlst;
    int ***lxyz;
    double **nfact;
    double *dfact;
public:
    AuxFunctions(int maxl)
    {
        int mlp1=maxl+1;
        nlst=new int[mlp1];
        for (int i=0;i<=maxl;i++) {
            nlst[i]=(i+1)*(i+2)/2;
        }
        maxlst=nlst[maxl];
        lxyz=new int**[mlp1];
        for (int i=0;i<=maxl;i++) {
            int nls=nlst[i];
            lxyz[i]=new int*[nls];
            for (int lx=0;lx<nls;++lx) lxyz[i][lx]=new int[3];
            int k=0;
            for (int lx=i;lx>=0;lx--) {
                for (int ly=i-lx;ly>=0;ly--) {
                    lxyz[i][k][0]=lx;
                    lxyz[i][k][1]=ly;
                    lxyz[i][k][2]=i-lx-ly;
                    k++;
                }
            }
        }
        double dfact_x[]={   1.0, 1.0, 3.0, 15.0, 105.0, 945.0 };
        dfact=new double[5];
        memcpy(dfact,dfact_x,sizeof(double)*5);
        nfact=new double*[mlp1];
        for (int i=0;i<=maxl;i++) {
            int nls=nlst[i];
            nfact[i]=new double[nls];
            for (int j=0;j<nls;j++) {
                double dx=dfact[lxyz[i][j][0]];
                double dy=dfact[lxyz[i][j][1]];
                double dz=dfact[lxyz[i][j][2]];
                nfact[i][j]=1.0/sqrt(dx*dy*dz);
            }
        }
    };
    int maxLstates() const { return maxlst;};
    // number of L vectors for lvalue lv //
    int number_of_lstates(int lv) const {
        return nlst[lv];
    };
    // the L vector for lvalue lv and l state ls //
    const int* lvector(int lv,int ls) const {
        return lxyz[lv][ls];
    };
    // tensor for storing the L vectors for lvalue lv and l state ls //
    int Lxyz(int lv,int ls,int i) const {
        return lxyz[lv][ls][i];
    };
    // vector of L dependent normalization factors //
    const double* norm_factor_vec(int lv) const {
        return nfact[lv];
    };
    // L vector dependent normalization factor //
    double normalization_factor(int lv,int ls) const {
        return nfact[lv][ls];
    };
    double angfactor(int lx) const { return dfact[lx];};
};
#endif
