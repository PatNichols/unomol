#ifndef UNOMOL_TWOELEC_HPP
#define UNOMOL_TWOELEC_HPP
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include "Util.hpp"
#include "Basis.hpp"
#include "Rys.hpp"
#include "Stopwatch.hpp"
#include "cache.hpp"
#include "MD_Dfunction.hpp"
#include "MD_Rfunction.hpp"
#define MAXFILESIZE 1073741824UL

namespace unomol {

#define UNO_MASK 0xF
#define UNO_SHIFT 4U
#define UNO_SHIFT2 8U

struct TwoInts {
    double val;
    int i,j,k,l;
};

struct ShellQuartet {
    double ab2,cd2;
    const double *a;
    const double *b;
    const double *c;
    const double *d;
    int npr1,lv1;
    const double *al1;
    const double *co1;
    int npr2,lv2;
    const double *al2;
    const double *co2;
    int npr3,lv3;
    const double *al3;
    const double *co3;
    int npr4,lv4;
    const double *al4;
    const double *co4;
    double ab[3];
    double cd[3];
    double * norms;
    unsigned int *lstates;
    unsigned int maxints;
    int len;
    bool sw12,sw34;
     
    ShellQuartet(int maxl) {
        int maxlst = ((maxl+1)*(maxl+2))/2;
        maxints = maxlst * maxlst * maxlst * maxlst;
        norms = new double[maxints];
        lstates = new unsigned int[maxints];
        len = 0;
    }

    ~ShellQuartet() {
        delete [] lstates;
        delete [] norms;
    }

    void assign_one(const Shell& sh,const double *rc) noexcept
    {
        npr1=sh.number_of_prims();
        lv1=sh.Lvalue();
        al1=sh.alf_ptr();
        co1=sh.cof_ptr();
        a=rc;
    }
    void assign_two(const Shell& sh,const double *rc) noexcept
    {
        npr2=sh.number_of_prims();
        lv2=sh.Lvalue();
        al2=sh.alf_ptr();
        co2=sh.cof_ptr();
        b=rc; 
        ab[0] = a[0] - b[0];
        ab[1] = a[1] - b[1];
        ab[2] = a[2] - b[2];
        ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2]*ab[2];
    }

    void assign_three(const Shell& sh,const double *rc) noexcept
    {
        npr3=sh.number_of_prims();
        lv3=sh.Lvalue();
        al3=sh.alf_ptr();
        co3=sh.cof_ptr();
        c=rc;
    }
    void assign_four(const Shell& sh,const double *rc) noexcept
    {
        npr4=sh.number_of_prims();
        lv4=sh.Lvalue();
        al4=sh.alf_ptr();
        co4=sh.cof_ptr();
        d=rc; 
        cd[0] = c[0] - d[0];
        cd[1] = c[1] - d[1];
        cd[2] = c[2] - d[2];
        cd2 = cd[0] * cd[0] + cd[1] * cd[1] + cd[2]*cd[2];
    }

    void swap_12() noexcept
    {
      sw12 = false;
      if ( lv1 < lv2) {
            sw12 = true;
                int it=npr1;
                npr1=npr2;
                npr2=it;
                it=lv1;
                lv1=lv2;
                lv2=it;
                const double * dp=al1;
                al1=al2;
                al2=dp;
                dp=co1;
                co1=co2;
                co2=dp;
                dp=a;
                a=b;
                b=dp;            
           
            ab[0] = -ab[0];
            ab[1] = -ab[1];
            ab[2] = -ab[2];
      }
    }

    void swap_34() noexcept
    {
      sw34 = false;
      if ( lv3 < lv4) {
            sw34 = true;
                int it=npr3;
                npr3=npr4;
                npr4=it;
                it=lv3;
                lv3=lv4;
                lv4=it;
                const double * dp=al3;
                al3=al4;
                al4=dp;
                dp=co3;
                co3=co4;
                co4=dp;
                dp=c;
                c=d;
                d=dp;            
           
            cd[0] = -cd[0];
            cd[1] = -cd[1];
            cd[2] = -cd[2];
      }
    }
 
    void unswap_12() noexcept
    {
      if (sw12) {
          sw12 = false;
                npr1=npr2;
                lv1=lv2;
                al1=al2;
                co1=co2;
                a=b;          
      }
    }

    void unswap_34() noexcept
    {
      if (sw34) {
          sw34 = false;
          npr3 = npr4;
          lv3 = lv4;
          al3 = al4;
          co3 = co4;
          c = d;
      }
    } 

   int precalculate(const int *offs,
            const AuxFunctions& aux,
            TwoInts * sints)
    {
/*
                    int nls1 = nls[0];
                    int nls2 = nls[1];
                    int nls3 = nls[2];
                    int nls4 = nls[3];
*/
                    int ir0 = offs[0];
                    int jr0 = offs[1];
                    int kr0 = offs[2];
                    int lr0 = offs[3];
                    int lv1_,lv2_,lv3_,lv4_;
                    if ( sw12 ) {
                        lv1_ = lv2;
                        lv2_ = lv1;
                    }else{
                        lv1_ = lv1;
                        lv2_ = lv2;
                    }
                    if ( sw34) {
                        lv3_ = lv4;
                        lv4_ = lv3;
                    }else{
                        lv3_ = lv3;
                        lv4_ = lv4;
                    }
                    const int nls1 = aux.number_of_lstates(lv1_);
                    const int nls2 = aux.number_of_lstates(lv2_);
                    const int nls3 = aux.number_of_lstates(lv3_);
                    const int nls4 = aux.number_of_lstates(lv4_);
                    int knt = 0;
                    for (int ils=0; ils<nls1;++ils) {
                        int ir = ir0 + ils;
                        for (int jls=0; jls<nls2;++jls) {
                            int jr = jr0 + jls;
                            if ( jr > ir ) break;
                            for (int kls=0; kls<nls3;++kls) {
                                int kr = kr0 + kls;
                                if ( kr > ir) break;
                                for (int lls=0; lls<nls4;++lls) {
                                    int lr = lr0 + lls;
                                    if ( lr > kr || ( ir == kr && lr > jr) ) break;
                                    (sints+knt)->val=0.0;
                                    (sints+knt)->i=(unsigned int)ir;
                                    (sints+knt)->j=(unsigned int)jr;
                                    (sints+knt)->k=(unsigned int)kr;
                                    (sints+knt)->l=(unsigned int)lr;
                                    unsigned int l12 = (ils<<UNO_SHIFT) + jls;
                                    if (sw12) l12=(jls<<UNO_SHIFT)+ils;
                                    unsigned int l34=(kls<<UNO_SHIFT)+lls;
                                    if (sw34) l34=(lls<<UNO_SHIFT)+kls;
                                    lstates[knt]=(l12<<UNO_SHIFT2)+l34;
                                    norms[knt] =
                                         aux.normalization_factor(lv1_,ils)*
                                         aux.normalization_factor(lv2_,jls)*
                                         aux.normalization_factor(lv3_,kls)*
                                         aux.normalization_factor(lv4_,lls);                     
                                    ++knt;
                                }
                            }
                        }
                    }
                    len = knt;
                    return len;
      }

};


struct MDInts {
    MD_Dfunction dx12;
    MD_Dfunction dy12;
    MD_Dfunction dz12;
    MD_Dfunction dx34;
    MD_Dfunction dy34;
    MD_Dfunction dz34;
    MD_Rfunction rfun;

    MDInts(int maxl):dx12(maxl),dy12(maxl),dz12(maxl),dx34(maxl),dy34(maxl),dz34(maxl),rfun(maxl) {
        std::cerr << "MD ints\n";
    }
};

class TwoElectronInts {
  public:
    TwoElectronInts() = delete;

    TwoElectronInts(const Basis& basis,
                                 int start_shell,const string& base_str):
        start(start_shell),rank(0),psize(1),cache(base_str,8192L*1048576L)
    {
        calculate(basis);
    }

    ~TwoElectronInts() {
    }

    void calculate(const Basis& base);

    void recalculate(const Basis& base) {
        calculate(base);
    }

    void formGmatrix(const double* Pmat,double *Gmat);

    void formGmatrix(const double* PmatA,const double* PmatB,
                     double* GmatA,double* GmatB);
  private:
    putils::Cache cache;
    int start;
    int rank,psize;
};


}
#endif
