#include <cstdlib>
#include <string>
#include "Basis.hpp"
#include "TwoElectronInts.hpp"
#include "RHF.hpp"
#include "UHF.hpp"

int main(int argc,char **argv)
    {
    string label("MINTS.DAT");
    unomol::Basis bas;
    unomol::TwoElectronInts t(bas,0,label);
    int nelec=bas.number_of_electrons();
    if (nelec%2)
        {
        unomol::UnRestrictedHartreeFock uhf(&bas,&t);
        uhf.findEnergy();
        if (bas.int_flags(1)) uhf.FiniteFieldAnalysis();
        if (bas.int_flags(0)) uhf.findPolarizationPotential();
        return EXIT_SUCCESS;
        }
    else
        {
        unomol::RestrictedHartreeFock rhf(&bas,&t);
        rhf.findEnergy();
        if (bas.int_flags(1)) rhf.FiniteFieldAnalysis();
        if (bas.int_flags(0)) rhf.findPolarizationPotential();
        return EXIT_SUCCESS;
        }
    }
