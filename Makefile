COMPILER=g++
FLAGS=-O3 -Wall -pedantic -fopenmp -std=c++11

MPI_COMPILER=mpixlC
MPI_FLAGS=-O3 -Wall -pedantic -fopenmp -std=c++11

all: main main-mpi

main: main.cpp
	$(COMPILER) $(FLAGS) main.cpp -o main

main-mpi: main_mpi.cpp
	$(MPI_COMPILER) $(MPI_FLAGS) main_mpi.cpp -o main_mpi

submit-polus: main-mpi
	rm -rf output.txt error.txt
	mpisubmit.pl -p 5 --stdout output.txt --stderr error.txt ./main_mpi
