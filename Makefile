CFLAGS=-g -Ofast -fPIC -ffast-math -mtune=native -Wno-unused-result 
CXXFLAGS=-g -Ofast -fPIC -ffast-math -mtune=native -Wno-unused-result 

SER_FILES= Util.o SymmPack.o Rys.o Moments.o OneElectronInts.o FField.o TwoElectronInts.o GDPMInts.o Unomol.o

MPI_FILES= Util.o SymmPack.o Rys.o Moments.o OneElectronInts.o FField.o TwoElectronIntsMPI.o GDPMInts.o UnomolMPI.o

all: 	$(SER_FILES)
	$(CXX) $(CXXFLAGS) -o Unomol $(SER_FILES)
	$(CXX) -o input input.cc
	$(CXX) -o test/srder test/srder.cc

mpilib:
	rm -f *.o
	$(MPICXX) $(CXXFLAGS) -c Util.cpp
	$(MPICXX) $(CXXFLAGS) -c SymmPack.cpp
	$(MPICXX) $(CXXFLAGS) -c Rys.cpp
	$(MPICXX) $(CXXFLAGS) -c OneElectronInts.cpp
	$(MPICXX) $(CXXFLAGS) -c FField.cpp
	$(MPICXX) $(CXXFLAGS) -c Moments.cpp
	$(MPICXX) $(CXXFLAGS) -c GDPMInts.cpp
	$(MPICXX) $(CXXFLAGS) -DUNOMOL_MPI_ABI -c TwoElectronIntsMPI.cpp
	$(MPICXX) $(CXXFLAGS) -c UnomolMPI.cc
	$(MPICXX) $(CXXFLAGS) -o UnomolMPI $(MPI_FILES)
		
clean:
	rm -fr *.o Unomol UnomolMPI input test/srder
	
