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
    const double *al1;
    const double *co1;
    const double *al2;
    const double *co2;
    const double *al3;
    const double *co3;
    const double *al4;
    const double *co4;
    double * norms;
    unsigned int *lstates;
    int npr1,lv1;
    int npr2,lv2;
    int npr3,lv3;
    int npr4,lv4;
    int len;
    int switch12,switch34;
    
    ShellQuartet(int maxl) {
        int maxlst = ((maxl+1)*(maxl+2))/2;
        int maxints = maxlst * maxlst * maxlst * maxlst;
        norms = new double[maxints];
        lstates = new unsigned int[maxints];
        switch12 = 0;
        switch34 = 0;
        len = 0;
    }

    ~ShellQuartet() {
        delete [] lstates;
        delete [] norms;
    }

    inline void assign1(const Shell& sh, const Center * center) noexcept
    {
      npr1 = sh.number_of_prims();
      lv1 = sh.Lvalue();
      al1 = sh.alf_ptr();
      co1 = sh.cof_ptr();
      a = (center+sh.center())->r_vec();
    }
    inline void assign3(const Shell& sh, const Center * center) noexcept
    { 
      npr3 = sh.number_of_prims();
      lv3 = sh.Lvalue();
      al3 = sh.alf_ptr();
      co3 = sh.cof_ptr();
      c = (center+sh.center())->r_vec();
    }
    inline void assign2(const Shell& sh, const Center * center) noexcept
    {
      npr2 = sh.number_of_prims();
      lv2 = sh.Lvalue();
      al2 = sh.alf_ptr();
      co2 = sh.cof_ptr();
      b = (center+sh.center())->r_vec();
      ab2 = dist_sqr(a,b);
    }
    inline void assign4(const Shell& sh, const Center * center) noexcept
    { 
      npr4 = sh.number_of_prims();
      lv4 = sh.Lvalue();
      al4 = sh.alf_ptr();
      co4 = sh.cof_ptr();
      d = (center+sh.center())->r_vec();
      cd2 = dist_sqr(c,d);
    }
/*
    inline void assign2(const Shell& sh, const Center * center) noexcept
    {
      lv2 = sh.Lvalue();
      if ( lv2 <= lv1) {
        npr2 = sh.number_of_prims();
        al2 = sh.alf_ptr();
        co2 = sh.cof_ptr();
        b = (center+sh.center())->r_vec();
      }else{
        b = a;
        al2 = al1;
        co2 = co1;
        npr2 = npr1;
        lv2 = lv1;
        npr1 = sh.number_of_prims();
        lv1 = sh.Lvalue();
        int cn1 = sh.center();
        al1 = sh.alf_ptr();
        co1 = sh.cof_ptr();
        a = (center+cn1)->r_vec();
        switch12 = 1;
      }
      ab2 = dist_sqr(a,b);
    }
    inline void assign4(const Shell& sh, const Center * center) noexcept
    {
      lv4 = sh.Lvalue();
      if ( lv4 <= lv3) {
        npr4 = sh.number_of_prims();
        lv4 = sh.Lvalue();
        al4 = sh.alf_ptr();
        co4 = sh.cof_ptr();
        d = (center+sh.center())->r_vec();
      }else{
        d = c;
        al4 = al3;
        co4 = co3;
        npr4 = npr3;
        lv4 = lv3;
        npr1 = sh.number_of_prims();
        lv3 = sh.Lvalue();
        int cn3 = sh.center();
        al3 = sh.alf_ptr();
        co3 = sh.cof_ptr();
        c = (center+cn3)->r_vec();
        switch34 = 1;
      }
      cd2 = dist_sqr(c,d);
    }
*/
    inline void unswitch12() noexcept
    {
      if ( switch12) {
        al1 = al2;
        co1 = co2;
        a = b;
        npr1 = npr2;
        lv1 = lv2;
      }
      switch12 = 0;
    }
    inline void unswitch34() noexcept
    {
      if ( switch34) {
        al3 = al4;
        co3 = co4;
        c = d;
        npr3 = npr4;
        lv3 = lv4;
      }
      switch34 = 0;
    }
    inline int precalculate( 
      TwoInts * sints,
      const AuxFunctions& aux,
      int off1,
      int off2,
      int off3,
      int off4) noexcept
    {
      len = 0;
      const int nls1 = aux.number_of_lstates(lv1);
      const int nls2 = aux.number_of_lstates(lv2);
      const int nls3 = aux.number_of_lstates(lv3);
      const int nls4 = aux.number_of_lstates(lv4);
      for (unsigned int ils=0;ils<nls1;++ils)
      {
        unsigned int ir1 = off1 + ils;
        for (unsigned int jls=0;jls<nls2;++jls)
        {
          unsigned int ir2 = off2 + jls;
          if ( ir2 > ir1) break;
          for (unsigned int kls=0;kls<nls3;++kls)
          {
            unsigned int ir3 = off3 + kls;
            if ( ir3 > ir1) break;
            for (unsigned int lls=0;lls<nls4;++lls)
            {
              unsigned int ir4 = off4 + lls;
              if ( ir4 > ir3) break;
              if ( ir1 == ir3 && ir4 > ir2) break;
              sints[len].val = 0.0;
              sints[len].i = ir1;
              sints[len].j = ir2;
              sints[len].k = ir3;
              sints[len].l = ir4;
              lstates[len] = (ils << 12) + ( jls << 8) + (kls << 4) + lls;
              norms[len] = aux.normalization_factor(lv1,ils) *
                            aux.normalization_factor(lv2,jls) *
                            aux.normalization_factor(lv3,kls) *
                            aux.normalization_factor(lv4,lls);
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

    void directFormGMatrix(const double *Pmat, double *Gmat,
        const Basis& base);
  private:
    putils::Cache cache;
    int start;
    int rank,psize;
};


}
#endif
