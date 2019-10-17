#include <cstdlib>
#include <string>
#include "Basis.cc"
#include "TwoElectronInts.cc"
#include "RHF.cc"
#include "UHF.cc"

int main()
{
    string label("MINTS.DAT");
    Basis bas;
    TwoElectronInts t(bas,0,label);
    int nelec=bas.number_of_electrons();
    if (nelec%2) {
        UnRestrictedHartreeFock uhf(&bas,&t);
        uhf.findEnergy();
        if (bas.int_flags(1)) uhf.FiniteFieldAnalysis();
        if (bas.int_flags(0)) uhf.findPolarizationPotential();
        return EXIT_SUCCESS;
    }else{
        RestrictedHartreeFock rhf(&bas,&t);
        rhf.findEnergy();
        if (bas.int_flags(1)) rhf.FiniteFieldAnalysis();
        if (bas.int_flags(0)) rhf.findPolarizationPotential();
        return EXIT_SUCCESS;
    }
}
