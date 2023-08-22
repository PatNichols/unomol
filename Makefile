
CFLAGS= -Ofast -ffast-math -mtune=native -mavx2 --std=c99
CXXFLAGS= -Ofast -O3 -ffast-math -mtune=native --std=c++17 -mavx2

#CXXFLAGS=-g -O1 -fPIC -fsanitize=address --std=c++17
#CFLAGS=-g -O1 -fPIC -fsanitize=address


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
	
