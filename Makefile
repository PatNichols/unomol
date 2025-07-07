CFLAGS=-g -O2 -fPIC -ffast-math -mtune=native
CXXFLAGS=-g -O2 -fPIC -ffast-math -mtune=native --std=c++20  
#CXXFLAGS=-g -O -fPIC -ggdb --std=c++17  

SER_FILES= putils_c.o Util.o SymmPack.o Rys.o Moments.o OneElectronInts.o FField.o TwoElectronInts.o GDPMInts.o Unomol.o

MPI_FILES= putils_c.o Util.o SymmPack.o Rys.o Moments.o OneElectronInts.o FField.o TwoElectronIntsMPI.o GDPMInts.o UnomolMPI.o

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<


all: 	$(SER_FILES)
	$(CXX) $(CXXFLAGS) -o Unomol $(SER_FILES)
	$(CXX) -o input input.cc
	$(CXX) -o test/srder test/srder.cc

mpilib:
	rm -f *.o
	$(MPICC) $(CFLAGS) -c putils_c.c
	$(MPICXX) $(CXXFLAGS) -c Util.cpp
	$(MPICXX) $(CXXFLAGS) -c SymmPack.cpp
	$(MPICXX) $(CXXFLAGS) -c Rys.cpp
	$(MPICXX) $(CXXFLAGS) -c OneElectronInts.cpp
	$(MPICXX) $(CXXFLAGS) -c FField.cpp
	$(MPICXX) $(CXXFLAGS) -c Moments.cpp
	$(MPICXX) $(CXXFLAGS) -c GDPMInts.cpp
	$(MPICXX) $(CXXFLAGS) -c TwoElectronIntsMPI.cpp
	$(MPICXX) $(CXXFLAGS) -c UnomolMPI.cc
	$(MPICXX) $(CXXFLAGS) -o UnomolMPI $(MPI_FILES)
		
clean:
	rm -fr *.o Unomol UnomolMPI input test/srder
	
