#pragma once

#include <iostream>
#include <cassert>
#include "Entities.h"

const char X_AXIS = 'x';
const char Y_AXIS = 'y';
const char Z_AXIS = 'z';

struct Volume {
    int xmin;
    int xmax;
    int ymin;
    int ymax;
    int zmin;
    int zmax;

    int dx;
    int dy;
    int dz;
    int size;
};

std::ostream& operator<<(std::ostream& os, const Volume &volume) {
    os << "[" << volume.xmin << ", " << volume.xmax << "] x [";
    os << volume.ymin << ", " << volume.ymax << "] x [";
    os << volume.zmin << ", " << volume.zmax << "]";
    return os;
}

class MPIHyperbolicEquationSolver {
    VolumeSize L; // параметры параллелепипеда
    double T; // время
    BoundaryConditionTypes bt; // граничные условия
    int N; // размер пространственной сетки
    int K; // размер временной сетки
    double hx, hy, hz; // шаги пространственной сетки
    double tau; // шаг временной сетки
    int rank; // номер процесса
    int size; // количество процессов

    int layerSize; // размер слоя
    Volume volume; // рабочий параллелепипед
    std::map<int, Volume> sendNeighbours; // соседи (процесс -> область соседа)
    std::map<int, Volume> recvNeighbours; // соседи (процесс -> область соседа)

    Volume MakeVolume(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) const;
    std::vector<double> PackVolume(const std::vector<double> &u, Volume volume) const;
    double FindValue(const std::vector<double> &u, int i, int j, int k, const std::vector<std::vector<double>> &u_recv) const;

    void SplitGrid(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes);
    bool IsInside(int xmin1, int xmax1, int ymin1, int ymax1, int xmin2, int xmax2, int ymin2, int ymax2) const;
    bool GetNeighbours(Volume v1, Volume v2, Volume &neighbour) const;
    void FillNeighbours(const std::vector<Volume> &volumes);

    int LocalIndex(int i, int j, int k) const;
    double GetBoundaryValue(int i, int j, int k, double t) const;

    void FillBoundaryValues(std::vector<double> &u, double t) const; // заполнение граничными значениями
    void FillInitialValues(std::vector<double> &u0, std::vector<double> &u1) const; // заполнение начальных условий

    double LaplaceOperator(const std::vector<double> &u, int i, int j, int k, const std::vector<std::vector<double>> &u_recv) const; // оператор Лапласа
    double EvaluateError(const std::vector<double> &u, double t) const; // оценка погрешности на слое
public:
    MPIHyperbolicEquationSolver(VolumeSize L, double T, int N, int K, BoundaryConditionTypes bt, int rank, int size);

    double AnalyticalSolve(double x, double y, double z, double t) const; // аналитическое решение
    double Phi(double x, double y, double z) const; // начальные условия

    void Solve(int maxSteps = 20, const char *numericalPath = NULL, const char *analyticalPath = NULL); // решение
    void PrintParams() const; // вывод параметров
};

MPIHyperbolicEquationSolver::MPIHyperbolicEquationSolver(VolumeSize L, double T, int N, int K, BoundaryConditionTypes bt, int rank, int size) {
    this->L = L;
    this->T = T;

    this->N = N;
    this->K = K;

    this->hx = L.x / N;
    this->hy = L.y / N;
    this->hz = L.z / N;
    this->tau = T / K;

    this->bt = bt;
    this->layerSize = (N + 1) * (N + 1) * (N + 1);

    this->rank = rank;
    this->size = size;
}

Volume MPIHyperbolicEquationSolver::MakeVolume(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) const {
    int dx = xmax - xmin + 1;
    int dy = ymax - ymin + 1;
    int dz = zmax - zmin + 1;
    int size = dx * dy * dz;

    return { xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, size };
}

std::vector<double> MPIHyperbolicEquationSolver::PackVolume(const std::vector<double> &u, Volume v) const {
    std::vector<double> packed(v.size);

    for (int i = v.xmin; i <= v.xmax; i++) {
        for (int j = v.ymin; j <= v.ymax; j++) {
            for (int k = v.zmin; k <= v.zmax; k++) {
                int index = (i - v.xmin) * v.dy * v.dz + (j - v.ymin) * v.dz + (k - v.zmin);
                packed[index] = u[LocalIndex(i, j, k)];
            }
        }
    }

    return packed;
}

void MPIHyperbolicEquationSolver::SplitGrid(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes) {
    if (size == 1) {
        volumes.push_back(MakeVolume(xmin, xmax, ymin, ymax, zmin, zmax));
        return;
    }

    if (size % 2 == 1) {
        if (axis == X_AXIS) {
            int x = xmin + (xmax - xmin) / size;

            volumes.push_back(MakeVolume(xmin, x, ymin, ymax, zmin, zmax));
            xmin = x + 1;
            axis = Y_AXIS;
        }
        else if (axis == Y_AXIS) {
            int y = ymin + (ymax - ymin) / size;
            volumes.push_back(MakeVolume(xmin, xmax, ymin, y, zmin, zmax));
            ymin = y + 1;
            axis = Z_AXIS;
        }
        else {
            int z = zmin + (zmax - zmin) / size;
            volumes.push_back(MakeVolume(xmin, xmax, ymin, ymax, zmin, z));
            zmin = z + 1;
            axis = X_AXIS;
        }

        size--;
    }

    if (axis == X_AXIS) {
        int x = (xmin + xmax) / 2;
        SplitGrid(xmin, x, ymin, ymax, zmin, zmax, size / 2, Y_AXIS, volumes);
        SplitGrid(x + 1, xmax, ymin, ymax, zmin, zmax, size / 2, Y_AXIS, volumes);
    }
    else if (axis == Y_AXIS) {
        int y = (ymin + ymax) / 2;
        SplitGrid(xmin, xmax, ymin, y, zmin, zmax, size / 2, Z_AXIS, volumes);
        SplitGrid(xmin, xmax, y + 1, ymax, zmin, zmax, size / 2, Z_AXIS, volumes);
    }
    else {
        int z = (zmin + zmax) / 2;
        SplitGrid(xmin, xmax, ymin, ymax, zmin, z, size / 2, X_AXIS, volumes);
        SplitGrid(xmin, xmax, ymin, ymax, z + 1, zmax, size / 2, X_AXIS, volumes);
    }
}

bool MPIHyperbolicEquationSolver::IsInside(int xmin1, int xmax1, int ymin1, int ymax1, int xmin2, int xmax2, int ymin2, int ymax2) const {
    return xmin2 <= xmin1 && xmax1 <= xmax2 && ymin2 <= ymin1 && ymax1 <= ymax2;
}

bool MPIHyperbolicEquationSolver::GetNeighbours(Volume v1, Volume v2, Volume &neighbour) const {
    if (v1.xmin == v2.xmax + 1 || v2.xmin == v1.xmax + 1) {
        int x = v1.xmin == v2.xmax + 1 ? v1.xmin : v1.xmax;

        if (IsInside(v1.ymin, v1.ymax, v1.zmin, v1.zmax, v2.ymin, v2.ymax, v2.zmin, v2.zmax)) {
            neighbour = MakeVolume(x, x, v1.ymin, v1.ymax, v1.zmin, v1.zmax);
            return true;
        }

        if (IsInside(v2.ymin, v2.ymax, v2.zmin, v2.zmax, v1.ymin, v1.ymax, v1.zmin, v1.zmax)) {
            neighbour = MakeVolume(x, x, v2.ymin, v2.ymax, v2.zmin, v2.zmax);
            return true;
        }

        return false;
    }

    if (v1.ymin == v2.ymax + 1 || v2.ymin == v1.ymax + 1) {
        int y = v1.ymin == v2.ymax + 1 ? v1.ymin : v1.ymax;

        if (IsInside(v1.xmin, v1.xmax, v1.zmin, v1.zmax, v2.xmin, v2.xmax, v2.zmin, v2.zmax)) {
            neighbour = MakeVolume(v1.xmin, v1.xmax, y, y, v1.zmin, v1.zmax);
            return true;
        }

        if (IsInside(v2.xmin, v2.xmax, v2.zmin, v2.zmax, v1.xmin, v1.xmax, v1.zmin, v1.zmax)) {
            neighbour = MakeVolume(v2.xmin, v2.xmax, y, y, v2.zmin, v2.zmax);
            return true;
        }

        return false;
    }

    if (v1.zmin == v2.zmax + 1 || v2.zmin == v1.zmax + 1) {
        int z = v1.zmin == v2.zmax + 1 ? v1.zmin : v1.zmax;

        if (IsInside(v1.xmin, v1.xmax, v1.ymin, v1.ymax, v2.xmin, v2.xmax, v2.ymin, v2.ymax)) {
            neighbour = MakeVolume(v1.xmin, v1.xmax, v1.ymin, v1.ymax, z, z);
            return true;
        }

        if (IsInside(v2.xmin, v2.xmax, v2.ymin, v2.ymax, v1.xmin, v1.xmax, v1.ymin, v1.ymax)) {
            neighbour = MakeVolume(v2.xmin, v2.xmax, v2.ymin, v2.ymax, z, z);
            return true;
        }

        return false;
    }

    return false;
}

void MPIHyperbolicEquationSolver::FillNeighbours(const std::vector<Volume> &volumes) {
    for (int i = 0; i < size; i++) {
        Volume neighbour;

        if (i == rank)
            continue;

        if (GetNeighbours(volume, volumes[i], neighbour)) {
            sendNeighbours[i] = neighbour;
        }

        if (GetNeighbours(volumes[i], volume, neighbour)) {
            recvNeighbours[i] = neighbour;
        }
    }
}

// аналитическое решение
double MPIHyperbolicEquationSolver::AnalyticalSolve(double x, double y, double z, double t) const {
    double at = M_PI * sqrt(1 / (L.x * L.x) + 1 / (L.y * L.y) + 1 / (L.z * L.z));

    return sin(M_PI / L.x * x) * sin(2 * M_PI / L.y * y) * sin(3 * M_PI / L.z * z) * cos(at * t + 2*M_PI);
}

// начальные условия
double MPIHyperbolicEquationSolver::Phi(double x, double y, double z) const {
    return AnalyticalSolve(x, y, z, 0);
}

int MPIHyperbolicEquationSolver::LocalIndex(int i, int j, int k) const {
    return (i - volume.xmin) * volume.dy * volume.dz + (j - volume.ymin) * volume.dz + (k - volume.zmin);
}

// TODO: numerical boundary values
double MPIHyperbolicEquationSolver::GetBoundaryValue(int i, int j, int k, double t) const {
    if (i == 0 || i == N) {
        if (bt.x == BoundaryConditionType::FirstKind)
            return 0;

        return AnalyticalSolve(i * hx, j * hy, k * hz, t);
    }

    if (j == 0 || j == N) {
        if (bt.y == BoundaryConditionType::FirstKind)
            return 0;

        return AnalyticalSolve(i * hx, j * hy, k * hz, t);
    }

    if (k == 0 || k == N) {
        if (bt.z == BoundaryConditionType::FirstKind)
            return 0;

        return AnalyticalSolve(i * hx, j * hy, k * hz, t);
    }

    throw "not boundary value";
}

// заполнение граничными значениями
void MPIHyperbolicEquationSolver::FillBoundaryValues(std::vector<double> &u, double t) const {
    if (volume.xmin == 0) {
        #pragma omp parallel for collapse(2)
        for (int i = volume.ymin; i <= volume.ymax; i++)
            for (int j = volume.zmin; j <= volume.zmax; j++)
                u[LocalIndex(volume.xmin, i, j)] = GetBoundaryValue(volume.xmin, i, j, t);
    }

    if (volume.xmax == N) {
        #pragma omp parallel for collapse(2)
        for (int i = volume.ymin; i <= volume.ymax; i++)
            for (int j = volume.zmin; j <= volume.zmax; j++)
                u[LocalIndex(volume.xmax, i, j)] = GetBoundaryValue(volume.xmax, i, j, t);
    }

    if (volume.ymin == 0) {
        #pragma omp parallel for collapse(2)
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.zmin; j <= volume.zmax; j++)
                u[LocalIndex(i, volume.ymin, j)] = GetBoundaryValue(i, volume.ymin, j, t);
    }

    if (volume.ymax == N) {
        #pragma omp parallel for collapse(2)
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.zmin; j <= volume.zmax; j++)
                u[LocalIndex(i, volume.ymax, j)] = GetBoundaryValue(i, volume.ymax, j, t);
    }

    if (volume.zmin == 0) {
        #pragma omp parallel for collapse(2)
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.ymin; j <= volume.ymax; j++)
                u[LocalIndex(i, j, volume.zmin)] = GetBoundaryValue(i, j, volume.zmin, t);
    }

    if (volume.zmax == N) {
        #pragma omp parallel for collapse(2)
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.ymin; j <= volume.ymax; j++)
                u[LocalIndex(i, j, volume.zmax)] = GetBoundaryValue(i, j, volume.zmax, t);
    }
}

void MPIHyperbolicEquationSolver::FillInitialValues(std::vector<double> &u0, std::vector<double> &u1) const {
    FillBoundaryValues(u0, 0);
    FillBoundaryValues(u1, tau);

    int xmin = std::max(volume.xmin, 1);
    int xmax = std::min(volume.xmax, N - 1);

    int ymin = std::max(volume.ymin, 1);
    int ymax = std::min(volume.ymax, N - 1);

    int zmin = std::max(volume.zmin, 1);
    int zmax = std::min(volume.zmax, N - 1);

    #pragma omp parallel for collapse(3)
    for (int i = xmin; i <= xmax; i++)
        for (int j = ymin; j <= ymax; j++)
            for (int k = zmin; k <= zmax; k++)
                u0[LocalIndex(i, j, k)] = Phi(i * hx, j * hy, k * hz);

    std::vector<std::vector<double>> u0_recv;

    for (auto it = sendNeighbours.begin(); it != sendNeighbours.end(); it++) {
        std::vector<double> packed = PackVolume(u0, it->second);
        MPI_Send(packed.data(), packed.size(), MPI_DOUBLE, it->first, 0, MPI_COMM_WORLD);
    }

    for (auto it = recvNeighbours.begin(); it != recvNeighbours.end(); it++) {
        std::vector<double> packed(it->second.size);
        MPI_Recv(packed.data(), packed.size(), MPI_DOUBLE, it->first, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        u0_recv.push_back(packed);
    }

    #pragma omp parallel for collapse(3)
    for (int i = xmin; i <= xmax; i++) {
        for (int j = ymin; j <= ymax; j++) {
            for (int k = zmin; k <= zmax; k++) {
                u1[LocalIndex(i, j, k)] = u0[LocalIndex(i, j, k)] + tau * tau / 2 * LaplaceOperator(u0, i, j, k, u0_recv);
            }
        }
    }
}

// TODO: refactore
double MPIHyperbolicEquationSolver::FindValue(const std::vector<double> &u, int i, int j, int k, const std::vector<std::vector<double>> &u_recv) const {
    int ind = 0;

    for (auto it = recvNeighbours.begin(); it != recvNeighbours.end(); it++) {
        Volume v = it->second;

        if (i < v.xmin || i > v.xmax || j < v.ymin || j > v.ymax || k < v.zmin || k > v.zmax) {
            ind++;
            continue;
        }

        int index = (i - v.xmin) * v.dy * v.dz + (j - v.ymin) * v.dz + (k - v.zmin);
        return u_recv[ind][index];

    }

    return u[LocalIndex(i, j, k)];
}

// оператор Лапласа
double MPIHyperbolicEquationSolver::LaplaceOperator(const std::vector<double> &u, int i, int j, int k, const std::vector<std::vector<double>> &u_recv) const {
    double dx = (FindValue(u, i - 1, j, k, u_recv) - 2 * u[LocalIndex(i, j, k)] + FindValue(u, i + 1, j, k, u_recv)) / (hx * hx);
    double dy = (FindValue(u, i, j - 1, k, u_recv) - 2 * u[LocalIndex(i, j, k)] + FindValue(u, i, j + 1, k, u_recv)) / (hy * hy);
    double dz = (FindValue(u, i, j, k - 1, u_recv) - 2 * u[LocalIndex(i, j, k)] + FindValue(u, i, j, k + 1, u_recv)) / (hz * hz);

    return dx + dy + dz;
}

// оценка погрешности на слое
double MPIHyperbolicEquationSolver::EvaluateError(const std::vector<double> &u, double t) const {
    double localError = 0;

    #pragma omp parallel for collapse(3) reduction(max: localError)
    for (int i = volume.xmin; i <= volume.xmax; i++)
        for (int j = volume.ymin; j <= volume.ymax; j++)
            for (int k = volume.zmin; k <= volume.zmax; k++)
                localError = std::max(localError, fabs(u[LocalIndex(i, j, k)] - AnalyticalSolve(i * hx, j * hy, k * hz, t)));

    double error;
    MPI_Reduce(&localError, &error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return error;
}

// решение
void MPIHyperbolicEquationSolver::Solve(int maxSteps, const char *numericalPath, const char *analyticalPath) {
    std::vector<Volume> volumes;
    SplitGrid(0, N, 0, N, 0, N, size, X_AXIS, volumes);
    volume = volumes[rank];

    FillNeighbours(volumes);

    std::vector<double> u0(volume.size);
    std::vector<double> u1(volume.size, -199);

    FillInitialValues(u0, u1);

    double error0 = EvaluateError(u0, 0);
    double error1 = EvaluateError(u1, tau);

    if (rank == 0) {
        std::cout << "Layer 0 max error: " << std::setprecision(15) << error0 << std::endl;
        std::cout << "Layer 1 max error: " << std::setprecision(15) << error1 << std::endl;
    }
}

// вывод параметров
void MPIHyperbolicEquationSolver::PrintParams() const {
    std::cout << "Processors: " << size << std::endl << std::endl;

    std::cout << "Lx: " << L.x << std::endl;
    std::cout << "Ly: " << L.y << std::endl;
    std::cout << "Lz: " << L.z << std::endl;
    std::cout << "T: " << T << std::endl << std::endl;

    std::cout << "N: " << N << std::endl;
    std::cout << "K: " << K << std::endl << std::endl;

    std::cout << "x_i = i*" << hx << ", i = 0..." << N << std::endl;
    std::cout << "y_i = i*" << hy << ", i = 0..." << N << std::endl;
    std::cout << "z_i = i*" << hz << ", i = 0..." << N << std::endl;
    std::cout << "t_i = i*" << tau << ", i = 0..." << K << std::endl << std::endl;

    std::cout << "Type of boundary condition along axis X: " << bt.x << std::endl;
    std::cout << "Type of boundary condition along axis Y: " << bt.y << std::endl;
    std::cout << "Type of boundary condition along axis Z: " << bt.z << std::endl;
}
