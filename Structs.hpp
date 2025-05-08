#ifndef UNOMOL_STRUCTS_HPP
#define UNOMOL_STRUCTS_HPP

namespace unomol {

struct ShellPairData {
    double ab2;
    const double *a,*b;
    int npr1,lv1;
    const double *al1,*co1;
    int npr2,lv2;
    const double *al2,*co2;
    unsigned int *lstates;
    double * norms;
    int * orbs;
    unsigned int len;
    int maxlst;

    ShellPairData() = delete;
        
    ShellPairData(int maxl)
    {
        if ( maxl < 2 ) maxl = 2;
        maxlst = ((maxl+1)*(maxl+2)) / 2;
        maxlst = maxlst * maxlst;
//        std::cerr << "sp maxlst = " << maxlst;
        lstates = new unsigned int[maxlst];
        norms = new double[maxlst];
        orbs = new int[maxlst]; 
    }

    ~ShellPairData()
    {
        delete [] orbs;
        delete [] norms;
        delete [] lstates;
    }

    void assign_one( const Shell & sh, const double *rc) noexcept {
        npr1 = sh.number_of_prims();
        lv1 = sh.Lvalue();
        al1 = sh.alf_ptr();
        co1 = sh.cof_ptr();
        a = rc;
    }
    void assign_two( const Shell & sh, const double *rc) noexcept {
        npr2 = sh.number_of_prims();
        lv2 = sh.Lvalue();
        al2 = sh.alf_ptr();
        co2 = sh.cof_ptr();
        b = rc;
        ab2 = 0.0;
        for (int i=0;i<3;++i) {
            double dab = a[i] - b[i];
            ab2 += dab * dab;
        }
    }

    int precalculate(int off1,int off2,const AuxFunctions& aux, double *f)
    {
        len = 0;
        int nls1 = aux.number_of_lstates(lv1);
        int nls2 = aux.number_of_lstates(lv2);
        for (int ls1=0;ls1<nls1;++ls1) {
            int ir = off1 + ls1;
            int iir = ir * ( ir +  1 ) / 2;
            for (int ls2=0;ls2<nls2;++ls2) {
                int jr = off2 + ls2;
                if ( jr > ir ) break;
                if ( len >= maxlst ) {
                    std::cerr << " len too large!\n";
                    exit(-1);
                }
                norms[len] = aux.normalization_factor(lv1,ls1) * aux.normalization_factor(lv2,ls2);
                lstates[len] = (ls1 << 4) + ls2;
                orbs[len] = iir + jr;
                if ( f ) {
                    if ( ir != jr ) {
                        f[len] = 1.0;
                    } else {
                        f[len] = 2.0;
                    }
                }
                ++len;
            }
        }
        return len;
    }

    int precalculate(int off1,int off2,const AuxFunctions& aux)
    {
        len = 0;
        int nls1 = aux.number_of_lstates(lv1);
        int nls2 = aux.number_of_lstates(lv2);
        for (int ls1=0;ls1<nls1;++ls1) {
            int ir = off1 + ls1;
            int iir = ir * ( ir +  1 ) / 2;
            for (int ls2=0;ls2<nls2;++ls2) {
                int jr = off2 + ls2;
                if ( jr > ir ) break;
                if ( len >= maxlst ) {
                    std::cerr << " len too large!\n";
                    exit(-1);
                }
                norms[len] = aux.normalization_factor(lv1,ls1) * aux.normalization_factor(lv2,ls2);
                lstates[len] = (ls1 << 4) + ls2;
                orbs[len] = iir + jr;
                ++len;
            }
        }
        return len;
    }

};

struct MomInts {
    double dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz;
    int ijr;
};

struct Moments
{
    double dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz;
};

}
#endif

