CFLAGS=-g -O3 -fPIC -ffast-math -march=native
CXXFLAGS=-g -O3 -fPIC -ffast-math -march=native
SER_FILES = Unomol.cc 

MPI_FILES= UnomolMPI.cc 

all: 	
	$(MPICXX) $(CFLAGS) -o UnomolMPI UnomolMPI.cc
	$(CXX) $(CFLAGS) -o Unomol Unomol.cc
	$(CXX) -o input input.cc
	$(CXX) -o test/srder test/srder.cc
	
clean:
	rm -fr Unomol UnomolMPI input test/srder
	
