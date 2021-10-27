#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include <mpi.h>

#include "ArgumentParser.hpp"
#include "MPIHyperbolicEquationSolver.hpp"

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
    ArgumentParser parser;

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

    MPIHyperbolicEquationSolver solver(arguments.L, arguments.T, arguments.N, arguments.K, arguments.bt, rank, size);

    if (arguments.debug && rank == 0) {
        cout << "Readed parameters: " << endl;
        solver.PrintParams(arguments.outputPath);
        cout << endl;
    }

    solver.Solve(arguments.steps, arguments.outputPath, arguments.numericalPath, arguments.analyticalPath);

    MPI_Finalize();
}