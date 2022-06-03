#include <cstdlib>
#include <string>
#define UNOMOL_MPI_ABI
#include <mpi.h>
#include "Basis.hpp"
#include "TwoElectronInts.hpp"
#include "RHF_MPI.hpp"
#include "UHF_MPI.hpp"

int main(int argc,char **argv) {
    MPI_Init(&argc,&argv);
    int rank,nproc;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    char buff[128];
    sprintf(buff,"MINTS_%04d.DAT",rank);
    std::string label(buff);
    unomol::Basis bas;
    unomol::TwoElectronInts t(bas,0,label);
    int nelec=bas.number_of_electrons();
    if (nelec%2) {
        unomol::UnRestrictedHartreeFockMPI uhf(&bas,&t);
        uhf.findEnergy();
        if (bas.int_flags(1)) uhf.FiniteFieldAnalysis();
        if (bas.int_flags(0)) uhf.findPolarizationPotential();
        MPI_Finalize();
        return EXIT_SUCCESS;
    } else {
        unomol::RestrictedHartreeFockMPI rhf(&bas,&t);
        rhf.findEnergy();
        if (bas.int_flags(1)) rhf.FiniteFieldAnalysis();
        if (bas.int_flags(0)) rhf.findPolarizationPotential();
        MPI_Finalize();
        return EXIT_SUCCESS;
    }
}
