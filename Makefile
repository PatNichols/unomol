CFLAGS=-g -O3 -fPIC -ffast-math -mtune=native -Wno-unused-result
CXXFLAGS=-g -O3 -fPIC -ffast-math -mtune=native -Wno-unused-result

SER_FILES= Util.o SymmPack.o Rys.o Moments.o OneElectronInts.o FField.o TwoElectronInts.o GDPMInts.o Unomol.o

all: 	$(SER_FILES)
	$(CXX) $(CXXFLAGS) -o Unomol $(SER_FILES)
	$(CXX) -o input input.cc
	$(CXX) -o test/srder test/srder.cc
	
clean:
	rm -fr *.o Unomol UnomolMPI input test/srder
	
