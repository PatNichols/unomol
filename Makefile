CXX = g++
CFLAGS = -O3 -funroll-loops -ffast-math -march=native -mavx2 -fno-exceptions 

FILES =	AuxFunctions.cc \
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
	$(MPI_CXX) $(CFLAGS) -o UnomolMPI UnomolMPI.cc

all: $(FILES)
	$(CXX) $(CFLAGS) -o Unomol Unomol.cc
	g++ -o input input.cc
	g++ -o test/srder test/srder.cc
	
clean:
	rm -fr Unomol UnomolMPI
	
