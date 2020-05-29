CXX = g++
CFLAGS= -g -O3 -fPIC -ffast-math -march=native
CXXFLAGS= -g -O3 -fPIC -ffast-math -march=native
SER_FILES = AuxFunctions.cc \
	Basis.cc \
	Center.cc \
	FField.cc \
	GDPMInts.cc \
	MD_Dfunction.cc \
	MD_Rfunction.cc \
	Moments.cc \
	OneElectronInts.cc \
	pconfig.h \
	RHF.cc \
	Rys.cc \
	Shell.cc \
	Structs.cc \
	SymmPack.cc \
	TwoElectronInts.cc \
	Unomol.cc \
	UHF.cc \
	Util.cc 

MPI_FILES= AuxFunctions.cc \
	Basis.cc \
	Center.cc \
	FField.cc \
	GDPMInts.cc \
	MD_Dfunction.cc \
	MD_Rfunction.cc \
	Moments.cc \
	OneElectronInts.cc \
	pconfig.h \
	RHF_MPI.cc \
	Rys.cc \
	Shell.cc \
	Structs.cc \
	SymmPack.cc \
	TwoElectronInts.cc \
	UnomolMPI.cc \
	UHF_MPI.cc \
	Util.cc 


mpif:	$(MPI_FILES)
	$(MPICXX) $(CFLAGS) -o UnomolMPI UnomolMPI.cc

all: 	$(SER_FILES)
	$(CXX) $(CFLAGS) -o Unomol Unomol.cc
	g++ -o input input.cc
	g++ -o test/srder test/srder.cc
	
clean:
	rm -fr Unomol UnomolMPI
	
