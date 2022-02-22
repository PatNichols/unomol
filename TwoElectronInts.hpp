#ifndef UNOMOL_TWOELEC_HPP
#define UNOMOL_TWOELEC_HPP
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include "pconfig.h"
#include "Util.hpp"
#include "Basis.hpp"
#include "Rys.hpp"
#include "Stopwatch.hpp"

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
    unsigned int *lstates;
    unsigned int len;

    ShellQuartet(int maxl) {
        int maxlst = ((maxl+1)*(maxl+2))/2;
        int maxints = maxlst * maxlst * maxlst * maxlst;
        lstates = new unsigned int[maxints];
    }

    ~ShellQuartet() {
        delete [] lstates;
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
                    int start_shell,const string& base_str) {
        basename = base_str;
        start = start_shell;
        rank = 0;
        psize = 1;
        calculate(basis);
    }
    ~TwoElectronInts() {
    }

    void calculate(const Basis& base);

    void recalculate(const Basis& base) {
        calculate(base);
    }

    void formGmatrix(const double* Pmat,double *Gmat) const;

    void formGmatrix(const double* PmatA,const double* PmatB,
                     double* GmatA,double* GmatB) const;
  private:
    std::string basename;
    int start;
    int rank,psize;
    std::vector<int> numints;
    int nfiles;

    FILE * create_ints_file(int file_no) const {
        if (file_no > 99999) {
	    fprintf(stderr,"max number of integral files exceeded!\n");
	    exit(-1);
        }
        char fname[127];
        int e = sprintf(fname,"%.5u_%.5u_%s",file_no,rank,basename.c_str());
        FILE * fp = fopen(fname,"w");
        if (fp) return fp;
        fprintf(stderr,"could not open the file %s\n",fname);
        exit(-1);
    }
    FILE * open_ints_file(int file_no) const {
        char fname[127];
        int e = sprintf(fname,"%.5u_%.5u_%s",file_no,rank,basename.c_str());
        FILE * fp = fopen(fname,"r");
        if (fp) return fp;
        fprintf(stderr,"could not open the file %s\n",fname);
        exit(-1);
    }
};


}
#endif
