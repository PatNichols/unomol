cmake_minimum_required( VERSION 3.10 FATAL_ERROR )

project( Unomol VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED True)

option(ENABLE_MPI "Build MPI version" OFF)

add_executable( Unomol
 	Unomol.cc
	FField.cpp
	GDPMInts.cpp
	Moments.cpp
	OneElectronInts.cpp
	Rys.cpp
	SymmPack.cpp
	TwoElectronInts.cpp
	Util.cpp
	putils_c.c)

target_compile_options( Unomol PUBLIC -O3 -ffast-math -mtune=native)
	
if (ENABLE_MPI)

	find_package(openmpi FATAL)

	add_executable( UnomolMPI
 		UnomolMPI.cc
		FField.cpp
		GDPMInts.cpp
		Moments.cpp
		OneElectronInts.cpp
		Rys.cpp
		SymmPack.cpp
		TwoElectronIntsMPI.cpp
		Util.cpp
		putils_c.c)
	target_compile_options( UNOMOL -O3 -ffast-math -mtune-native)

endif()

add_executable( input input.cc)

add_subdirectory(test)

add_subdirectory(scat)


