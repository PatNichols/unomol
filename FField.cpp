
#include "FField.hpp"

namespace unomol {

void FiniteFieldMatrix( const Basis& basis, double *Ham,
                        const double *Efield) {
    const int no=basis.number_of_orbitals();
    const int no2=no*(no+1)/2;
    MomInts mint;
    FILE *infile=open_file("RMOM.DAT");
    for (int i=0; i<no2; ++i) {
        Fread(&mint,sizeof(MomInts),1,infile);
        int ijr=mint.ijr;
        Ham[ijr] -= (Efield[0]*mint.dx+
                     Efield[1]*mint.dy+
                     Efield[2]*mint.dz);
    }
    fclose(infile);
}

}
