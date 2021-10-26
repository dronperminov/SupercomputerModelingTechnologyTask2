#pragma once

#include <iostream>
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
    std::map<int, Volume> neighbours; // соседи (процесс -> область соседа)

    void SplitGrid(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes);
    bool IsInside(int xmin1, int xmax1, int ymin1, int ymax1, int xmin2, int xmax2, int ymin2, int ymax2) const;
    bool GetNeighbours(Volume v1, Volume v2, Volume &neighbour) const;
    void FillNeighbours(const std::vector<Volume> &volumes);
public:
    MPIHyperbolicEquationSolver(VolumeSize L, double T, int N, int K, BoundaryConditionTypes bt, int rank, int size);

    void Solve(int maxSteps = 20, const char *numericalPath = NULL, const char *analyticalPath = NULL); // решение
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

void MPIHyperbolicEquationSolver::SplitGrid(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes) {
    if (size == 1) {
        volumes.push_back({ xmin, xmax, ymin, ymax, zmin, zmax });
        return;
    }

    if (size % 2 == 1) {
        if (axis == X_AXIS) {
            int x = xmin + (xmax - xmin) / size;

            volumes.push_back({ xmin, x, ymin, ymax, zmin, zmax });
            xmin = x + 1;
            axis = Y_AXIS;
        }
        else if (axis == Y_AXIS) {
            int y = ymin + (ymax - ymin) / size;
            volumes.push_back({ xmin, xmax, ymin, y, zmin, zmax });
            ymin = y + 1;
            axis = Z_AXIS;
        }
        else {
            int z = zmin + (zmax - zmin) / size;
            volumes.push_back({ xmin, xmax, ymin, ymax, zmin, z });
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
            neighbour = { x, x, v1.ymin, v1.ymax, v1.zmin, v1.zmax };
            return true;
        }

        if (IsInside(v2.ymin, v2.ymax, v2.zmin, v2.zmax, v1.ymin, v1.ymax, v1.zmin, v1.zmax)) {
            neighbour = { x, x, v2.ymin, v2.ymax, v2.zmin, v2.zmax };
            return true;
        }

        return false;
    }

    if (v1.ymin == v2.ymax + 1 || v2.ymin == v1.ymax + 1) {
        int y = v1.ymin == v2.ymax + 1 ? v1.ymin : v1.ymax;

        if (IsInside(v1.xmin, v1.xmax, v1.zmin, v1.zmax, v2.xmin, v2.xmax, v2.zmin, v2.zmax)) {
            neighbour = { v1.xmin, v1.xmax, y, y, v1.zmin, v1.zmax };
            return true;
        }

        if (IsInside(v2.xmin, v2.xmax, v2.zmin, v2.zmax, v1.xmin, v1.xmax, v1.zmin, v1.zmax)) {
            neighbour = { v2.xmin, v2.xmax, y, y, v2.zmin, v2.zmax };
            return true;
        }

        return false;
    }

    if (v1.zmin == v2.zmax + 1 || v2.zmin == v1.zmax + 1) {
        int z = v1.zmin == v2.zmax + 1 ? v1.zmin : v1.zmax;

        if (IsInside(v1.xmin, v1.xmax, v1.ymin, v1.ymax, v2.xmin, v2.xmax, v2.ymin, v2.ymax)) {
            neighbour = { v1.xmin, v1.xmax, v1.ymin, v1.ymax, z, z };
            return true;
        }

        if (IsInside(v2.xmin, v2.xmax, v2.ymin, v2.ymax, v1.xmin, v1.xmax, v1.ymin, v1.ymax)) {
            neighbour = { v2.xmin, v2.xmax, v2.ymin, v2.ymax, z, z };
            return true;
        }

        return false;
    }

    return false;
}

void MPIHyperbolicEquationSolver::FillNeighbours(const std::vector<Volume> &volumes) {
    for (int i = 0; i < size; i++) {
        Volume neighbour;

        if (i != rank && GetNeighbours(volume, volumes[i], neighbour)) {
            neighbours[i] = neighbour;
        }
    }
}

// решение
void MPIHyperbolicEquationSolver::Solve(int maxSteps, const char *numericalPath, const char *analyticalPath) {
    std::vector<Volume> volumes;
    SplitGrid(0, N, 0, N, 0, N, size, X_AXIS, volumes);
    volume = volumes[rank];

    FillNeighbours(volumes);

    std::cout << "rank " << rank << ": " << volume << std::endl;

    for (auto it = neighbours.begin(); it != neighbours.end(); it++) {
        std::cout << "rank " << rank << ": " << it->first << " (" << it->second << std::endl;
    }

    std::cout << std::endl;
}
