COMPILER=g++
FLAGS=-O3 -Wall -pedantic -fopenmp

MPI_COMPILER=mpixlC
MPI_FLAGS=-O3 -Wall -pedantic -fopenmp -std=c++11

GPU_COMPILER=nvcc
GPU_FLAGS=-O3 -std=c++11 -I/opt/ibm/spectrum_mpi/include -L/opt/ibm/spectrum_mpi/lib -lmpiprofilesupport -lmpi_ibm

N=256
K=2000

all: main main-mpi

main: main.cpp
	$(COMPILER) $(FLAGS) main.cpp -o main

main-mpi: main_mpi.cpp
	$(MPI_COMPILER) $(MPI_FLAGS) main_mpi.cpp -o main_mpi

main-gpu: main_gpu.cu
	$(GPU_COMPILER) $(GPU_FLAGS) main_gpu.cu -o main_gpu

test-mpi: test_mpi.cpp
	$(MPI_COMPILER) $(MPI_FLAGS) test_mpi.cpp -o test_mpi

test-gpu: test_gpu.cu
	$(GPU_COMPILER) $(GPU_FLAGS) test_gpu.cu -o test_gpu

submit-polus: main-mpi
	mpisubmit.pl -p 16 --stdout stdout.txt --stderr error.txt ./main_mpi -- -d -o output.txt -N $(N) -K $(K)

submit-bluegene: main-mpi
	mpisubmit.bg -n 16 -m smp --stdout stdout.txt --stderr error.txt ./main_mpi -- -d -o output.txt -N $(N) -K $(K)

submit-polus-gpu: main-gpu
	bsub -n 16 -gpu "num=2" -R "span[ptile=2]" -o stdout.txt -e error.txt OMP_NUM_THREADS=1 mpiexec ./main_gpu -d -o output.txt -N $(N) -K $(K)

submit-polus-test: test-mpi
	for N in 128 256 512; do \
		for L in 1 3.14159265358979; do \
			for SPLIT in blocks tapes product; do \
				for p in 1 2 4 8 10 16 20 32 40 64; do \
					mpisubmit.pl -p $$p -w 00:30 --stdout /dev/null --stderr /dev/null ./test_mpi -- test_L$$L\_$$N\_$$SPLIT.txt $$N $$L $$SPLIT; \
				done \
			done \
		done \
	done

submit-bluegene-test: test-mpi
	for N in 128 256 512; do \
		for L in 1 3.14159265358979; do \
			for SPLIT in blocks tapes product; do \
				mpisubmit.bg -n 256 -w 00:10:00 -m smp --stdout /dev/null --stderr /dev/null ./test_mpi -- test_L$$L\_$$N\_$$SPLIT.txt $$N $$L $$SPLIT ; \
				mpisubmit.bg -n 128 -w 00:15:00 -m smp --stdout /dev/null --stderr /dev/null ./test_mpi -- test_L$$L\_$$N\_$$SPLIT.txt $$N $$L $$SPLIT ; \
				mpisubmit.bg -n 64 -w 00:30:00 -m smp --stdout /dev/null --stderr /dev/null ./test_mpi -- test_L$$L\_$$N\_$$SPLIT.txt $$N $$L $$SPLIT; \
				for p in 32 16 8 4 2 1; do \
					mpisubmit.bg -n $$p -w 01:59:00 -m smp --stdout /dev/null --stderr /dev/null ./test_mpi -- test_L$$L\_$$N\_$$SPLIT.txt $$N $$L $$SPLIT ; \
				done \
			done \
		done \
	done

submit-polus-gpu-test: test-gpu
	for N in 128 256 512; do \
		for L in 1; do \
			for SPLIT in blocks tapes product; do \
				for p in 1 2 4 8 16 32 64; do \
					bsub -n $$p -W 00:30 -gpu "num=2" -J test_gpu_L$$L\_$$N\_$$SPLIT\_P$$p -R "span[ptile=2]" -o /dev/null -e /dev/null OMP_NUM_THREADS=1 mpiexec ./test_gpu test_gpu_L$$L\_$$N\_$$SPLIT.txt $$N $$L $$SPLIT; \
				done \
			done \
		done \
	done