COMPILER=g++
FLAGS=-O3 -Wall -pedantic -fopenmp -std=c++11

MPI_COMPILER=mpixlC
MPI_FLAGS=-O3 -Wall -pedantic -fopenmp -std=c++11

N=256
K=2000

all: main main-mpi

main: main.cpp
	$(COMPILER) $(FLAGS) main.cpp -o main

main-mpi: main_mpi.cpp
	$(MPI_COMPILER) $(MPI_FLAGS) main_mpi.cpp -o main_mpi

test-mpi: test_mpi.cpp
	$(MPI_COMPILER) $(MPI_FLAGS) test_mpi.cpp -o test_mpi

submit-polus: main-mpi
	mpisubmit.pl -p 16 --stdout stdout.txt --stderr error.txt ./main_mpi -- -d -o output.txt -N $(N) -K $(K)

submit-polus-test: test-mpi
	for N in 128 256 512; do \
		for p in 1 2 4 8 10 16 20 32 40 64; do \
			for L in 1 3.14159265358979; do \
				mpisubmit.pl -p $$p -w 00:20 --stdout /dev/null --stderr /dev/null ./test_mpi -- test_L$$L\_$$N.txt $$N $$L ; \
			done \
		done \
	done
