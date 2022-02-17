#ifndef UNOMOL_TWOELEC_HPP
#define UNOMOL_TWOELEC_HPP
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "pconfig.h"
#include "Util.hpp"
#include "Basis.hpp"
#include "Rys.hpp"
#ifdef _MPI_API_
#include <mpi.h>
#endif

#define MAXFILESIZE 134217728

namespace unomol
{

struct Sints
    {
    double val;
    int i,j,k,l;
    };

struct ShellQuartet
    {
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

    ShellQuartet(int maxl)
        {
        int maxlst = ((maxl+1)*(maxl+2))/2;
        int maxints = maxlst * maxlst * maxlst * maxlst;
        lstates = new unsigned int[maxints];
        }

    ~ShellQuartet()
        {
        delete [] lstates;
        }
    };

std::uint64_t find_max_files(int no_new,int no_old,int nprocs);

class TwoElectronInts
    {
    public:
        TwoElectronInts() = delete;
        TwoElectronInts(const Basis& basis,
                        int start_shell,const string& base_str)
            {
            strcpy(basename,base_str.c_str());
            start = start_shell;
            rank = 0;
            psize = 1;
#ifdef UNOMOL_MPI_ABI_
            MPI_Comm_rank(MPI_COMM_WORLD,&rank);
            MPI_Comm_size(MPI_COMM_WORLD,&psize);
#endif
            int maxfiles = find_max_files(basis.number_of_orbitals(),0,psize);
            numints = new int[maxfiles];
            calculate(basis);
            }
        ~TwoElectronInts() {
            delete [] numints;
        }

        void calculate(const Basis& base);

        void recalculate(const Basis& base)
            {
            calculate(base);
            }

        void formGmatrix(const double* Pmat,double *Gmat) const;

        void formGmatrix(const double* PmatA,const double* PmatB,
                         double* GmatA,double* GmatB) const;
    private:
        char basename[72];
        int start;
        int rank,psize;
        int *numints;
        int nfiles;
        void calc_two_electron_ints(const ShellQuartet& sq,
                                    const AuxFunctions& aux,
                                    Rys& rys,
                                    Sints* sints);
    };


}
#endif
