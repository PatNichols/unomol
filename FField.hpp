#ifndef UNOMOL_FFIELD_hpp_
#define UNOMOL_FFIELD_hpp_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "Basis.hpp"
#include "AuxFunctions.hpp"
#include "MD_Dfunction.hpp"
#include "MD_Rfunction.hpp"
#include "Util.hpp"
#include "Structs.hpp"
using namespace std;

namespace unomol
{

inline void FiniteFieldMatrix( const Basis& basis, double *Ham,
                               const double *Efield)
    {
    const int no=basis.number_of_orbitals();
    const int no2=no*(no+1)/2;
    MomInts mint;
    FILE *infile=open_file("RMOM.DAT");
    for (int i=0; i<no2; ++i)
        {
        fread(&mint,sizeof(MomInts),1,infile);
        int ijr=mint.ijr;
        Ham[ijr] -= (Efield[0]*mint.dx+
                     Efield[1]*mint.dy+
                     Efield[2]*mint.dz);
        }
    fclose(infile);
    }

}
#endif

