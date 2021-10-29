#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
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

    std::vector<Volume> sendNeighbours; // соседи на передачу
    std::vector<Volume> recvNeighbours; // соседи на приём
    std::vector<int> processNeighbours; // соседи процессы

    Volume MakeVolume(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) const;
    std::vector<double> PackVolume(const std::vector<double> &u, Volume volume) const;
    std::vector<std::vector<double>> SendRecvValues(const std::vector<double> &u) const; // отправка/получение соседних значений
    std::vector<double> SendRecvTotal(const std::vector<double> &u, const std::vector<Volume> &volumes) const; // отправка/получение общих начений

    void SplitGrid(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes);
    bool IsInside(int xmin1, int xmax1, int ymin1, int ymax1, int xmin2, int xmax2, int ymin2, int ymax2) const;
    bool GetNeighbours(Volume v1, Volume v2, Volume &neighbour) const;
    void FillNeighbours(const std::vector<Volume> &volumes);

    int Index(int i, int j, int k, Volume v) const;
    int LocalIndex(int i, int j, int k) const;
    double GetBoundaryValue(int i, int j, int k, double t) const;

    void FillBoundaryValues(std::vector<double> &u, double t) const; // заполнение граничными значениями
    void FillInitialValues(std::vector<double> &u0, std::vector<double> &u1) const; // заполнение начальных условий
    void FillNextLayer(const std::vector<double> &u0, const std::vector<double> &u1, std::vector<double> &u, double t);
    void FillAnalyticalValues(std::vector<double> &u, double t) const; // заполнение аналитических значений
    void FillDifferenceValues(std::vector<double> &u, double t) const; // заполнение разности значений

    double FindValue(const std::vector<double> &u, int i, int j, int k, const std::vector<std::vector<double>> &u_recv) const;
    double LaplaceOperator(const std::vector<double> &u, int i, int j, int k, const std::vector<std::vector<double>> &u_recv) const; // оператор Лапласа
    double EvaluateError(const std::vector<double> &u, double t) const; // оценка погрешности на слое
    void SaveValues(const std::vector<double> u, double t, const std::vector<Volume> &volumes, const char *filename) const; // сохранение слоя
public:
    MPIHyperbolicEquationSolver(VolumeSize L, double T, int N, int K, BoundaryConditionTypes bt, int rank, int size);

    double AnalyticalSolve(double x, double y, double z, double t) const; // аналитическое решение
    double Phi(double x, double y, double z) const; // начальные условия

    double Solve(SolveParams params, bool isSilent = false); // решение
    void PrintParams(const char *outputPath) const; // вывод параметров
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

// аналитическое решение
double MPIHyperbolicEquationSolver::AnalyticalSolve(double x, double y, double z, double t) const {
    double at = M_PI * sqrt(1 / (L.x * L.x) + 1 / (L.y * L.y) + 1 / (L.z * L.z));

    return sin(M_PI / L.x * x) * sin(2 * M_PI / L.y * y) * sin(3 * M_PI / L.z * z) * cos(at * t + 2*M_PI);
}

// начальные условия
double MPIHyperbolicEquationSolver::Phi(double x, double y, double z) const {
    return AnalyticalSolve(x, y, z, 0);
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

    #pragma omp parallel for collapse(3)
    for (int i = v.xmin; i <= v.xmax; i++)
        for (int j = v.ymin; j <= v.ymax; j++)
            for (int k = v.zmin; k <= v.zmax; k++)
                packed[Index(i, j, k, v)] = u[LocalIndex(i, j, k)];

    return packed;
}

// отправка/получение соседних значений
std::vector<std::vector<double>> MPIHyperbolicEquationSolver::SendRecvValues(const std::vector<double> &u) const {
    std::vector<std::vector<double>> u_recv(processNeighbours.size());

    for (auto i = 0; i < processNeighbours.size(); i++) {
        std::vector<double> packed = PackVolume(u, sendNeighbours[i]);
        u_recv[i] = std::vector<double>(recvNeighbours[i].size);

        std::vector<MPI_Request> requests(2);
        std::vector<MPI_Status> statuses(2);

        MPI_Isend(packed.data(), sendNeighbours[i].size, MPI_DOUBLE, processNeighbours[i], 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(u_recv[i].data(), recvNeighbours[i].size, MPI_DOUBLE, processNeighbours[i], 0, MPI_COMM_WORLD, &requests[1]);
        MPI_Waitall(2, requests.data(), statuses.data());
    }

    return u_recv;
}

// отправка/получение общих начений
std::vector<double> MPIHyperbolicEquationSolver::SendRecvTotal(const std::vector<double> &u, const std::vector<Volume> &volumes) const {
    if (rank != 0) {
        MPI_Request request;
        MPI_Status status;

        MPI_Isend(u.data(), volume.size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Waitall(1, &request, &status);
        return u;
    }

    std::vector<double> u_all(layerSize);
    Volume volumeAll = MakeVolume(0, N, 0, N, 0, N);

    for (int index = 0; index < size; index++) {
        std::vector<double> ui(volumes[index].size);

        if (index == rank) {
            ui = u;
        }
        else {
            std::vector<MPI_Request> requests(1);
            std::vector<MPI_Status> statuses(1);

            MPI_Irecv(ui.data(), volumes[index].size, MPI_DOUBLE, index, 0, MPI_COMM_WORLD, &requests[0]);
            MPI_Waitall(1, requests.data(), statuses.data());
        }

        for (int i = volumes[index].xmin; i <= volumes[index].xmax; i++) {
            for (int j = volumes[index].ymin; j <= volumes[index].ymax; j++) {
                for (int k = volumes[index].zmin; k <= volumes[index].zmax; k++) {
                    u_all[Index(i, j, k, volumeAll)] = ui[Index(i, j, k, volumes[index])];
                }
            }
        }
    }

    return u_all;
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
    sendNeighbours.clear();
    recvNeighbours.clear();
    processNeighbours.clear();

    for (int i = 0; i < size; i++) {
        if (i == rank)
            continue;

        Volume sendNeighbour;
        Volume recvNeighbour;

        if (!GetNeighbours(volume, volumes[i], sendNeighbour))
            continue;

        GetNeighbours(volumes[i], volume, recvNeighbour);
        processNeighbours.push_back(i);
        sendNeighbours.push_back(sendNeighbour);
        recvNeighbours.push_back(recvNeighbour);
    }
}

int MPIHyperbolicEquationSolver::Index(int i, int j, int k, Volume v) const {
    return (i - v.xmin) * v.dy * v.dz + (j - v.ymin) * v.dz + (k - v.zmin);
}

int MPIHyperbolicEquationSolver::LocalIndex(int i, int j, int k) const {
    return Index(i, j, k, volume);
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
    for (int i = xmin; i <= xmax; i++) {
        for (int j = ymin; j <= ymax; j++) {
            for (int k = zmin; k <= zmax; k++) {
                u0[LocalIndex(i, j, k)] = Phi(i * hx, j * hy, k * hz);
            }
        }
    }

    std::vector<std::vector<double>> u0_recv = SendRecvValues(u0);

    #pragma omp parallel for collapse(3)
    for (int i = xmin; i <= xmax; i++) {
        for (int j = ymin; j <= ymax; j++) {
            for (int k = zmin; k <= zmax; k++) {
                u1[LocalIndex(i, j, k)] = u0[LocalIndex(i, j, k)] + tau * tau / 2 * LaplaceOperator(u0, i, j, k, u0_recv);
            }
        }
    }
}

void MPIHyperbolicEquationSolver::FillNextLayer(const std::vector<double> &u0, const std::vector<double> &u1, std::vector<double> &u, double t) {
    int xmin = std::max(volume.xmin, 1);
    int xmax = std::min(volume.xmax, N - 1);

    int ymin = std::max(volume.ymin, 1);
    int ymax = std::min(volume.ymax, N - 1);

    int zmin = std::max(volume.zmin, 1);
    int zmax = std::min(volume.zmax, N - 1);

    std::vector<std::vector<double>> u_recv = SendRecvValues(u1);

    #pragma omp parallel for collapse(3)
    for (int i = xmin; i <= xmax; i++)
        for (int j = ymin; j <= ymax; j++)
            for (int k = zmin; k <= zmax; k++)
                u[LocalIndex(i, j, k)] = 2 * u1[LocalIndex(i, j, k)] - u0[LocalIndex(i, j, k)] + tau * tau * LaplaceOperator(u1, i, j, k, u_recv);

    FillBoundaryValues(u, t);
}

// заполнение аналитических значений
void MPIHyperbolicEquationSolver::FillAnalyticalValues(std::vector<double> &u, double t) const {
    #pragma omp parallel for collapse(3)
    for (int i = volume.xmin; i <= volume.xmax; i++)
        for (int j = volume.ymin; j <= volume.ymax; j++)
            for (int k = volume.zmin; k <= volume.zmax; k++)
                u[LocalIndex(i, j, k)] = AnalyticalSolve(i * hx, j * hy, k * hz, t);
}

// заполнение разности значений
void MPIHyperbolicEquationSolver::FillDifferenceValues(std::vector<double> &u, double t) const {
    #pragma omp parallel for collapse(3)
    for (int i = volume.xmin; i <= volume.xmax; i++) {
        for (int j = volume.ymin; j <= volume.ymax; j++) {
            for (int k = volume.zmin; k <= volume.zmax; k++) {
                u[LocalIndex(i, j, k)] = fabs(u[LocalIndex(i, j, k)] - AnalyticalSolve(i * hx, j * hy, k * hz, t));
            }
        }
    }
}

double MPIHyperbolicEquationSolver::FindValue(const std::vector<double> &u, int i, int j, int k, const std::vector<std::vector<double>> &u_recv) const {
    for (auto it = 0; it < processNeighbours.size(); it++) {
        Volume v = recvNeighbours[it];

        if (i < v.xmin || i > v.xmax || j < v.ymin || j > v.ymax || k < v.zmin || k > v.zmax)
            continue;

        return u_recv[it][Index(i, j, k, v)];
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

// сохранение слоя
void MPIHyperbolicEquationSolver::SaveValues(const std::vector<double> u, double t, const std::vector<Volume> &volumes, const char *filename) const {
    std::vector<double> u_all = SendRecvTotal(u, volumes);

    if (rank != 0)
        return;

    std::ofstream f(filename);

    f << "{" << std::endl;
    f << "    \"Lx\": " << L.x << ", " << std::endl;
    f << "    \"Ly\": " << L.y << ", " << std::endl;
    f << "    \"Lz\": " << L.z << ", " << std::endl;
    f << "    \"N\": " << N << ", " << std::endl;
    f << "    \"t\": " << t << ", " << std::endl;
    f << "    \"u\": [" << std::endl;

    bool wasPrinted = false;

    for (int i = 0; i < layerSize; i++) {
        if (wasPrinted) {
            f << ", " << std::endl;
        }
        else {
            wasPrinted = true;
        }

        f << "    " << u_all[i];
    }

    f << std::endl;
    f << "    ]" << std::endl;
    f << "}" << std::endl;

    f.close();
}

// решение
double MPIHyperbolicEquationSolver::Solve(SolveParams params, bool isSilent) {
    std::vector<Volume> volumes;
    SplitGrid(0, N, 0, N, 0, N, size, X_AXIS, volumes);
    volume = volumes[rank];

    FillNeighbours(volumes);

    std::vector<std::vector<double>> u(3, std::vector<double>(volume.size));

    FillInitialValues(u[0], u[1]);

    double error0 = EvaluateError(u[0], 0);
    double error1 = EvaluateError(u[1], tau);

    if (rank == 0 && !isSilent) {
        std::ofstream fout(params.outputPath ? params.outputPath : "output.txt", std::ios::app);
        fout << "Layer 0 max error: " << error0 << std::endl;
        fout << "Layer 1 max error: " << error1 << std::endl;
        fout.close();
    }

    for (int step = 2; step <= params.steps; step++) {
        FillNextLayer(u[(step + 1) % 3], u[(step + 2) % 3], u[step % 3], step * tau);

        double error = EvaluateError(u[step % 3], step * tau);

        if (rank == 0 && !isSilent) {
            std::ofstream fout(params.outputPath ? params.outputPath : "output.txt", std::ios::app);
            fout << "Layer " << step << " max error: " << error << std::endl;
            fout.close();
        }
    }

    if (params.numericalPath) {
        SaveValues(u[params.steps % 3], params.steps * tau, volumes, params.numericalPath);
    }

    if (params.differencePath) {
        FillDifferenceValues(u[params.steps % 3], params.steps * tau);
        SaveValues(u[params.steps % 3], params.steps * tau, volumes, params.differencePath);
    }

    if (params.analyticalPath) {
        FillAnalyticalValues(u[0], params.steps * tau);
        SaveValues(u[0], params.steps * tau, volumes, params.analyticalPath);
    }

    return EvaluateError(u[params.steps % 3], params.steps * tau);
}

// вывод параметров
void MPIHyperbolicEquationSolver::PrintParams(const char *outputPath) const {
    std::ofstream fout(outputPath ? outputPath : "output.txt", std::ios::app);

    fout << "Processors: " << size << std::endl << std::endl;

    fout << "Lx: " << L.x << std::endl;
    fout << "Ly: " << L.y << std::endl;
    fout << "Lz: " << L.z << std::endl;
    fout << "T: " << T << std::endl << std::endl;

    fout << "N: " << N << std::endl;
    fout << "K: " << K << std::endl << std::endl;

    fout << "x_i = i*" << hx << ", i = 0..." << N << std::endl;
    fout << "y_i = i*" << hy << ", i = 0..." << N << std::endl;
    fout << "z_i = i*" << hz << ", i = 0..." << N << std::endl;
    fout << "t_i = i*" << tau << ", i = 0..." << K << std::endl << std::endl;

    fout << "Type of boundary condition along axis X: " << bt.x << std::endl;
    fout << "Type of boundary condition along axis Y: " << bt.y << std::endl;
    fout << "Type of boundary condition along axis Z: " << bt.z << std::endl;

    fout.close();
}
