#include <iostream>
#include <fstream>
#include <iomanip>
#include <mpi.h>

#include "include/ArgumentParser.hpp"
#include "include/CudaMPIHyperbolicEquationSolver.hpp"

using namespace std;

void InitMPI(int &rank, int &size) {
    if (MPI_Init(NULL, NULL) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, -1);
        throw "MPI_Init failed";
    }

    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, -1);
        throw "MPI_Comm_size failed";
    }

    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        MPI_Abort(MPI_COMM_WORLD, -1);
        throw "MPI_Comm_rank failed";
    }
}

int main(int argc, char **argv) {
    ArgumentParser parser(true);

    if (argc == 2 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
        parser.Help();
        return 0;
    }

    Arguments arguments;
    int rank, size;

    try {
        arguments = parser.Parse(argc, argv);
        InitMPI(rank, size);
    }
    catch (const char *error) {
        return -1;
    }

    double t0 = MPI_Wtime();

    if (arguments.debug && rank == 0) {
        PrintParams(arguments.L, arguments.T, arguments.N, arguments.K, arguments.bt, arguments.split, size, arguments.solveParams.outputPath);
    }

    CudaMPIHyperbolicEquationSolve(arguments.L, arguments.T, arguments.N, arguments.K, arguments.bt, arguments.split, rank, size, arguments.solveParams, false);

    double t1 = MPI_Wtime();
    double delta = t1 - t0;

    double maxTime, minTime, avgTime;

    MPI_Reduce(&delta, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&delta, &minTime, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&delta, &avgTime, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::ofstream fout(arguments.solveParams.outputPath ? arguments.solveParams.outputPath : "output.txt", std::ios::app);
        fout << "Max time: " << maxTime << std::endl;
        fout << "Min time: " << minTime << std::endl;
        fout << "Average time: " << avgTime / size << std::endl;
        fout.close();
    }

    MPI_Finalize();
}