#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <mpi.h>

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

bool IsExists(const char *path) {
    ifstream fin(path);

    if (!fin)
        return false;

    fin.close();
    return true;
}

SplitType GetSplit(const char *arg) {
    if (!strcmp(arg, "blocks") || !strcmp(arg, "b"))
        return SplitType::Blocks;

    if (!strcmp(arg, "tapes") || !strcmp(arg, "t"))
        return SplitType::Tapes;

    if (!strcmp(arg, "product") || !strcmp(arg, "p"))
        return SplitType::Product;

    throw "incorrect split strategy";
}

int main(int argc, char **argv) {
    int rank, size;

    try {
        InitMPI(rank, size);
    }
    catch (const char *error) {
        return -1;
    }

    int N = atoi(argv[2]);
    double L = atof(argv[3]);
    SplitType split = argc > 4 ? GetSplit(argv[4]) : SplitType::Blocks;

    double T = 1;
    int K = 2000;
    int loops = N <= 128 ? 20 : 5;

    BoundaryConditionTypes bt = {
        BoundaryConditionType::FirstKind,
        BoundaryConditionType::PeriodicAnalytical,
        BoundaryConditionType::FirstKind
    };

    SolveParams params;
    params.steps = 20;
    params.outputPath = NULL;
    params.numericalPath = NULL;
    params.analyticalPath = NULL;
    params.differencePath = NULL;

    if (!IsExists(argv[1])) {
        ofstream fout(argv[1]);
        fout << "### Lx = Ly = Lz = " << L << ", N = " << N << ", K = " << K << ", разбиение: " << split << endl << endl;
        fout << "| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |" << endl;
        fout << "|                     :-: |               :-: |       :-: |         :-: |" << endl;
        fout.close();
    }

    double t0 = MPI_Wtime();
    double error = 0;

    for (int loop = 0; loop < loops; loop++) {
        error += CudaMPIHyperbolicEquationSolve({ L, L, L }, T, N, K, bt, split, rank, size, params, true);
    }

    double t1 = MPI_Wtime();
    double delta = t1 - t0;
    double maxTime;

    MPI_Reduce(&delta, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        ofstream fout(argv[1], ios::app);
        fout << "| " << setw(23) << size;
        fout << " | " << setw(17) << (maxTime / loops);
        fout << " | " << "         ";
        fout << " | " << setw(11) << (error / loops);
        fout << " |" << endl;
        fout.close();
    }

    MPI_Finalize();
}