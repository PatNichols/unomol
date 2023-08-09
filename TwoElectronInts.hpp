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

    void assign_one(const Shell& sh,const double *rc)
    {
        npr1=sh.number_of_prims();
        lv1=sh.Lvalue();
        al1=sh.alf_ptr();
        co1=sh.cof_ptr();
        a=rc;
    }
    void assign_two(const Shell& sh,const double *rc)
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

    void assign_three(const Shell& sh,const double *rc)
    {
        npr3=sh.number_of_prims();
        lv3=sh.Lvalue();
        al3=sh.alf_ptr();
        co3=sh.cof_ptr();
        c=rc;
    }
    void assign_four(const Shell& sh,const double *rc)
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

    void swap_12()
    {
      sw12 = false;
 /*
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
 */
    }

    void swap_34()
    {
      sw34 = false;
/*
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
*/
    }
 
    void unswap_12()
    {
/*
      if (sw12) {
          sw12 = false;
                npr1=npr2;
                lv1=lv2;
                al1=al2;
                co1=co2;
                a=b;          
      }
*/      
    }

    void unswap_34()
    {
      if (sw34) {
/*
          sw34 = false;
          npr3 = npr4;
          lv3 = lv4;
          al3 = al4;
          co3 = co4;
          c = d;
*/
      }
    } 

    int precalculate(const AuxFunctions& aux,const int offs[],TwoInts *sints)
    {
         int lv1_,lv2_,lv3_,lv4_;
         int lsh12,lsh34;
         if ( sw12 ) {
              lv1_ = lv2;
              lv2_ = lv1; 
         } else {
              lv1_ = lv1;
              lv2_ = lv2;         
         }
         if ( sw34 ) {  
              lv3_ = lv4;
              lv4_ = lv3; 
         } else {
              lv3_ = lv3;
              lv4_ = lv4;         
         }
         int nls1 = aux.number_of_lstates(lv1_);
         int nls2 = aux.number_of_lstates(lv2_);
         int nls3 = aux.number_of_lstates(lv3_);
         int nls4 = aux.number_of_lstates(lv4_);
         len = 0;
         for (int ls1=0;ls1<nls1;++ls1) {
              double d1 = aux.normalization_factor(lv1_,ls1);
              int ir = offs[0] + ls1;
              for (int ls2 = 0; ls2 < nls2; ++ls2) {
                 double d12 = d1 * aux.normalization_factor(lv2_,ls2);
                 int jr = offs[1] + ls2;
                 if ( jr > ir) break;
                 if ( sw12) {
                    lsh12 = (ls2 << UNO_SHIFT) + ls1;  
                 }else{
                    lsh12 = (ls1 << UNO_SHIFT) + ls2;
                 }
                 for (int ls3=0;ls3<nls3;++ls3) {
                     double d123 = d12 * aux.normalization_factor(lv3_,ls3);
                     int kr = offs[2] + ls3;
                     if ( kr > ir ) break;
                     for (int ls4=0;ls4<nls4;++ls4) {
//                         norms[len] = d123 * aux.normalization_factor(lv4,ls4);
                         int lr = offs[3] + ls4;
                         if ( ( lr > kr ) || ( (ir == kr) && ( lr > jr) ) ) break;
                         if ( sw34 ) {
                            lsh34 = ( ls4 << UNO_SHIFT ) + ls3;
                         } else {
                            lsh34 = ( ls3 << UNO_SHIFT ) + ls4;
                         }
                         lstates[len] = (lsh12 << UNO_SHIFT2) + lsh34;
                         norms[len] = d123 * aux.normalization_factor(lv4_,ls4); 
                         sints[len].val = 0.0;
                         sints[len].i = ir;
                         sints[len].j = jr;
                         sints[len].k = kr;
                         sints[len].l = lr;
                         ++len;
                     } 
                 } 
              }
         }    
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
