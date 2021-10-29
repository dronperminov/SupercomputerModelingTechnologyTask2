#include <iostream>
#include <fstream>
#include <iomanip>
#include <mpi.h>

#include "include/MPIHyperbolicEquationSolver.hpp"

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

    double T = 1;
    int K = 2000;
    int steps = 20;
    int loops = 5;

    BoundaryConditionTypes bt = {
        BoundaryConditionType::FirstKind,
        BoundaryConditionType::PeriodicAnalytical,
        BoundaryConditionType::FirstKind
    };

    if (!IsExists(argv[1])) {
        ofstream fout(argv[1]);
        fout << "### Lx = Ly = Lz = " << L << ", N = " << N << ", K = " << K << endl << endl;
        fout << "| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |" << endl;
        fout << "|                     :-: |               :-: |       :-: |         :-: |" << endl;
        fout.close();
    }

    double t0 = MPI_Wtime();
    double error = 0;

    for (int loop = 0; loop < loops; loop++) {
        MPIHyperbolicEquationSolver solver({ L, L, L }, T, N, K, bt, rank, size);
        error += solver.Solve(steps, NULL, NULL, NULL, NULL, true);
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