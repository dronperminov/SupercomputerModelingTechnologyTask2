#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "Entities.h"
#include "GridSplitter.hpp"

const bool USE_GPU = true;
const int BLOCK_SIZE = 256; // рамер GPU блока

// аналитическое решение
__host__ __device__ double AnalyticalSolve(double x, double y, double z, double t, VolumeSize L) {
    double at = M_PI * sqrt(1 / (L.x * L.x) + 1 / (L.y * L.y) + 1 / (L.z * L.z));

    return sin(M_PI / L.x * x) * sin(2 * M_PI / L.y * y) * sin(3 * M_PI / L.z * z) * cos(at * t + 2*M_PI);
}

// начальные условия
__host__ __device__ double Phi(double x, double y, double z, VolumeSize L) {
    return AnalyticalSolve(x, y, z, 0, L);
}

__host__ __device__ int Index(int i, int j, int k, Volume v) {
    return (i - v.xmin) * v.dy * v.dz + (j - v.ymin) * v.dz + (k - v.zmin);
}

__host__ __device__ double FindValue(int size, const double *u, int i, int j, int k, const double *u_recv, const Volume *recvNeighbours, Volume volume) {
    if (volume.xmin <= i && i <= volume.xmax && volume.ymin <= j && j <= volume.ymax && volume.zmin <= k && k <= volume.zmax)
        return u[Index(i, j, k, volume)];

    int offset = 0;

    for (int index = 0; index < size; index++) {
        Volume v = recvNeighbours[index];

        if (i < v.xmin || i > v.xmax || j < v.ymin || j > v.ymax || k < v.zmin || k > v.zmax) {
            offset += recvNeighbours[index].size;
            continue;
        }

        return u_recv[offset + Index(i, j, k, v)];
    }

    return 42; // never
}

// оператор Лапласа
__host__ __device__ double LaplaceOperator(int size, const double *u, int i, int j, int k, const double *u_recv, const Volume *recvNeighbours, Volume volume, double hx, double hy, double hz) {
    double dx = (FindValue(size, u, i - 1, j, k, u_recv, recvNeighbours, volume) - 2 * u[Index(i, j, k, volume)] + FindValue(size, u, i + 1, j, k, u_recv, recvNeighbours, volume)) / (hx * hx);
    double dy = (FindValue(size, u, i, j - 1, k, u_recv, recvNeighbours, volume) - 2 * u[Index(i, j, k, volume)] + FindValue(size, u, i, j + 1, k, u_recv, recvNeighbours, volume)) / (hy * hy);
    double dz = (FindValue(size, u, i, j, k - 1, u_recv, recvNeighbours, volume) - 2 * u[Index(i, j, k, volume)] + FindValue(size, u, i, j, k + 1, u_recv, recvNeighbours, volume)) / (hz * hz);

    return dx + dy + dz;
}

__global__ void FillZeroLayerKernel(double *u0, int total, int xmin, int ymin, int zmin, int dy, int dz, Volume volume, double hx, double hy, double hz, VolumeSize L) {
    int index = blockIdx.x*blockDim.x + threadIdx.x;

    if (index >= total)
        return;

    int i = xmin + index / (dy*dz);
    int j = ymin + index % (dy*dz) / dz;
    int k = zmin + index % dz;

    u0[Index(i, j, k, volume)] = Phi(i * hx, j * hy, k * hz, L);
}

__global__ void EvaluateErrorKernel(double *u, double hx, double hy, double hz, double t, VolumeSize L, Volume volume, double *error) {
    int idx = threadIdx.x;
    double max = 0;

    for (int index = idx; index < volume.size; index += BLOCK_SIZE) {
        int i = volume.xmin + index / (volume.dy*volume.dz);
        int j = volume.ymin + index % (volume.dy*volume.dz) / volume.dz;
        int k = volume.zmin + index % volume.dz;
        double delta = fabs(u[Index(i, j, k, volume)] - AnalyticalSolve(i * hx, j * hy, k * hz, t, L));

        if (delta > max)
            max = delta;
    }

    __shared__ double r[BLOCK_SIZE];

    r[idx] = max;
    __syncthreads();

    for (int size = BLOCK_SIZE / 2; size > 0; size /= 2) {
        if (idx < size) {
            r[idx] = r[idx] > r[idx+size] ? r[idx] : r[idx+size];
        }

        __syncthreads();
    }

    if (idx == 0) {
        *error = r[0];
    }
}

__global__ void FillFirstLayerKernel(double *u1, const double *u0, const double *u0_recv, const Volume *recvNeighbours, int total, int size, int xmin, int ymin, int zmin, int dy, int dz, Volume volume, double hx, double hy, double hz, double tau) {
    int index = blockIdx.x*blockDim.x + threadIdx.x;

    if (index >= total)
        return;

    int i = xmin + index / (dy*dz);
    int j = ymin + index % (dy*dz) / dz;
    int k = zmin + index % dz;
    
    u1[Index(i, j, k, volume)] = u0[Index(i, j, k, volume)] + tau * tau / 2 * LaplaceOperator(size, u0, i, j, k, u0_recv, recvNeighbours, volume, hx, hy, hz);
}

__global__ void FillNextLayerKernel(double *u, const double *u0, const double *u1, const double *u_recv, const Volume *recvNeighbours, int total, int size, int xmin, int ymin, int zmin, int dy, int dz, Volume volume, double hx, double hy, double hz, double tau) {
    int index = blockIdx.x*blockDim.x + threadIdx.x;

    if (index >= total)
        return;

    int i = xmin + index / (dy*dz);
    int j = ymin + index % (dy*dz) / dz;
    int k = zmin + index % dz;

    u[Index(i, j, k, volume)] = 2 * u1[Index(i, j, k, volume)] - u0[Index(i, j, k, volume)] + tau * tau * LaplaceOperator(size, u1, i, j, k, u_recv, recvNeighbours, volume, hx, hy, hz);
}

void checkCUDAError(const char *msg) {
    cudaError_t err = cudaGetLastError();
    if (cudaSuccess != err) {
        std::cerr << "Cuda error: " << msg << ": " << cudaGetErrorString( err) << "." << std::endl;
        exit(-1);
    }                      
}

std::vector<double> PackVolume(const std::vector<double> &u, Volume v, Volume volume) {
    std::vector<double> packed(v.size);

    for (int i = v.xmin; i <= v.xmax; i++)
        for (int j = v.ymin; j <= v.ymax; j++)
            for (int k = v.zmin; k <= v.zmax; k++)
                packed[Index(i, j, k, v)] = u[Index(i, j, k, volume)];

    return packed;
}

// отправка/получение соседних значений
std::vector<double> SendRecvValues(const std::vector<double> &u, std::vector<int> processNeighbours, std::vector<Volume> sendNeighbours, std::vector<Volume> recvNeighbours, Volume volume) {
    std::vector<double> u_recv;
    int offset = 0;

    for (auto i = 0; i < processNeighbours.size(); i++) {
        std::vector<double> packed = PackVolume(u, sendNeighbours[i], volume);
        u_recv.insert(u_recv.end(), recvNeighbours[i].size, 0);

        std::vector<MPI_Request> requests(2);
        std::vector<MPI_Status> statuses(2);

        MPI_Isend(packed.data(), sendNeighbours[i].size, MPI_DOUBLE, processNeighbours[i], 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(u_recv.data() + offset, recvNeighbours[i].size, MPI_DOUBLE, processNeighbours[i], 0, MPI_COMM_WORLD, &requests[1]);
        MPI_Waitall(2, requests.data(), statuses.data());
        offset += recvNeighbours[i].size;
    }

    return u_recv;
}

// отправка/получение общих начений
std::vector<double> SendRecvTotal(const std::vector<double> &u, const std::vector<Volume> &volumes, int rank, int size, Volume volume, int layerSize, int N) {
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

bool IsInside(int xmin1, int xmax1, int ymin1, int ymax1, int xmin2, int xmax2, int ymin2, int ymax2) {
    return xmin2 <= xmin1 && xmax1 <= xmax2 && ymin2 <= ymin1 && ymax1 <= ymax2;
}

bool GetNeighbours(Volume v1, Volume v2, Volume &neighbour) {
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

void FillNeighbours(const std::vector<Volume> &volumes, std::vector<int> &processNeighbours, std::vector<Volume> &sendNeighbours, std::vector<Volume> &recvNeighbours, int rank, int size, Volume volume) {
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

double GetBoundaryValue(int i, int j, int k, double t, BoundaryConditionTypes bt, double hx, double hy, double hz, VolumeSize L, int N) {
    if (i == 0 || i == N) {
        if (bt.x == BoundaryConditionType::FirstKind)
            return 0;

        return AnalyticalSolve(i * hx, j * hy, k * hz, t, L);
    }

    if (j == 0 || j == N) {
        if (bt.y == BoundaryConditionType::FirstKind)
            return 0;

        return AnalyticalSolve(i * hx, j * hy, k * hz, t, L);
    }

    if (k == 0 || k == N) {
        if (bt.z == BoundaryConditionType::FirstKind)
            return 0;

        return AnalyticalSolve(i * hx, j * hy, k * hz, t, L);
    }

    throw "not boundary value";
}

// заполнение граничными значениями
void FillBoundaryValues(std::vector<double> &u, double t, BoundaryConditionTypes bt, double hx, double hy, double hz, VolumeSize L, Volume volume, int N) {
    if (volume.xmin == 0) {
        for (int i = volume.ymin; i <= volume.ymax; i++)
            for (int j = volume.zmin; j <= volume.zmax; j++)
                u[Index(volume.xmin, i, j, volume)] = GetBoundaryValue(volume.xmin, i, j, t, bt, hx, hy, hz, L, N);
    }

    if (volume.xmax == N) {
        for (int i = volume.ymin; i <= volume.ymax; i++)
            for (int j = volume.zmin; j <= volume.zmax; j++)
                u[Index(volume.xmax, i, j, volume)] = GetBoundaryValue(volume.xmax, i, j, t, bt, hx, hy, hz, L, N);
    }

    if (volume.ymin == 0) {
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.zmin; j <= volume.zmax; j++)
                u[Index(i, volume.ymin, j, volume)] = GetBoundaryValue(i, volume.ymin, j, t, bt, hx, hy, hz, L, N);
    }

    if (volume.ymax == N) {
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.zmin; j <= volume.zmax; j++)
                u[Index(i, volume.ymax, j, volume)] = GetBoundaryValue(i, volume.ymax, j, t, bt, hx, hy, hz, L, N);
    }

    if (volume.zmin == 0) {
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.ymin; j <= volume.ymax; j++)
                u[Index(i, j, volume.zmin, volume)] = GetBoundaryValue(i, j, volume.zmin, t, bt, hx, hy, hz, L, N);
    }

    if (volume.zmax == N) {
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.ymin; j <= volume.ymax; j++)
                u[Index(i, j, volume.zmax, volume)] = GetBoundaryValue(i, j, volume.zmax, t, bt, hx, hy, hz, L, N);
    }
}

void FillInitialValues(std::vector<double> &u0, std::vector<double> &u1, BoundaryConditionTypes bt, double hx, double hy, double hz, double tau, VolumeSize L, int N, const std::vector<int> &processNeighbours, const std::vector<Volume> &sendNeighbours, const std::vector<Volume> &recvNeighbours, Volume volume) {
    FillBoundaryValues(u0, 0, bt, hx, hy, hz, L, volume, N);
    FillBoundaryValues(u1, tau, bt, hx, hy, hz, L, volume, N);

    int xmin = std::max(volume.xmin, 1);
    int xmax = std::min(volume.xmax, N - 1);

    int ymin = std::max(volume.ymin, 1);
    int ymax = std::min(volume.ymax, N - 1);

    int zmin = std::max(volume.zmin, 1);
    int zmax = std::min(volume.zmax, N - 1);

    int dx = xmax - xmin + 1;
    int dy = ymax - ymin + 1;
    int dz = zmax - zmin + 1;
    int total = dx*dy*dz;

    if (!USE_GPU) {
        for (int index = 0; index < total; index++) {
            int i = xmin + index / (dy*dz);
            int j = ymin + index % (dy*dz) / dz;
            int k = zmin + index % dz;

            u0[Index(i, j, k, volume)] = Phi(i * hx, j * hy, k * hz, L);
        }

        std::vector<double> u0_recv = SendRecvValues(u0, processNeighbours, sendNeighbours, recvNeighbours, volume);

        for (int index = 0; index < total; index++) {
            int i = xmin + index / (dy*dz);
            int j = ymin + index % (dy*dz) / dz;
            int k = zmin + index % dz;
            
            u1[Index(i, j, k, volume)] = u0[Index(i, j, k, volume)] + tau * tau / 2 * LaplaceOperator(recvNeighbours.size(), u0.data(), i, j, k, u0_recv.data(), recvNeighbours.data(), volume, hx, hy, hz);
        }
    }
    else {
        int nblocks = (total + BLOCK_SIZE - 1) / BLOCK_SIZE;

        double *u0Device;
        cudaMalloc((void**) &u0Device, u0.size() * sizeof(double));
        FillZeroLayerKernel<<<nblocks, BLOCK_SIZE>>>(u0Device, total, xmin, ymin, zmin, dy, dz, volume, hx, hy, hz, L);
        cudaMemcpy(u0.data(), u0Device, u0.size() * sizeof(double), cudaMemcpyDeviceToHost);

        std::vector<double> u0_recv = SendRecvValues(u0, processNeighbours, sendNeighbours, recvNeighbours, volume);

        double *u0recvDevice;
        double *u1Device;
        Volume *recvNeighboursDevice;

        cudaMalloc((void**) &u1Device, u1.size() * sizeof(double));
        cudaMalloc((void**) &u0recvDevice, u0_recv.size() * sizeof(double));
        cudaMalloc((void**) &recvNeighboursDevice, recvNeighbours.size() * sizeof(Volume));

        cudaMemcpy(u0recvDevice, u0_recv.data(), u0_recv.size() * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(recvNeighboursDevice, recvNeighbours.data(), recvNeighbours.size() * sizeof(Volume), cudaMemcpyHostToDevice);

        FillFirstLayerKernel<<<nblocks, BLOCK_SIZE>>>(u1Device, u0Device, u0recvDevice, recvNeighboursDevice, total, recvNeighbours.size(), xmin, ymin, zmin, dy, dz, volume, hx, hy, hz, tau);
        cudaMemcpy(u1.data(), u1Device, u1.size() * sizeof(double), cudaMemcpyDeviceToHost);

        cudaFree(u0Device);
        cudaFree(u0recvDevice);
        cudaFree(u1Device);
        cudaFree(recvNeighboursDevice);
    }
}

void FillNextLayer(const std::vector<double> &u0, const std::vector<double> &u1, std::vector<double> &u, double t, double hx, double hy, double hz, double tau, VolumeSize L, int N, BoundaryConditionTypes bt, const std::vector<int> &processNeighbours, const std::vector<Volume> &sendNeighbours, const std::vector<Volume> &recvNeighbours, Volume volume) {
    int xmin = std::max(volume.xmin, 1);
    int xmax = std::min(volume.xmax, N - 1);

    int ymin = std::max(volume.ymin, 1);
    int ymax = std::min(volume.ymax, N - 1);

    int zmin = std::max(volume.zmin, 1);
    int zmax = std::min(volume.zmax, N - 1);

    int dx = xmax - xmin + 1;
    int dy = ymax - ymin + 1;
    int dz = zmax - zmin + 1;
    int total = dx*dy*dz;

    std::vector<double> u_recv = SendRecvValues(u1, processNeighbours, sendNeighbours, recvNeighbours, volume);

    if (!USE_GPU) {
        for (int index = 0; index < total; index++) {
            int i = xmin + index / (dy*dz);
            int j = ymin + index % (dy*dz) / dz;
            int k = zmin + index % dz;

            u[Index(i, j, k, volume)] = 2 * u1[Index(i, j, k, volume)] - u0[Index(i, j, k, volume)] + tau * tau * LaplaceOperator(recvNeighbours.size(), u1.data(), i, j, k, u_recv.data(), recvNeighbours.data(), volume, hx, hy, hz);
        }
    }
    else {
        int nblocks = (total + BLOCK_SIZE - 1) / BLOCK_SIZE;

        double *u0Device;
        double *u1Device;
        double *uDevice;
        double *urecvDevice;
        Volume *recvNeighboursDevice;

        cudaMalloc((void**) &u0Device, u0.size() * sizeof(double));
        cudaMalloc((void**) &u1Device, u1.size() * sizeof(double));
        cudaMalloc((void**) &uDevice, u.size() * sizeof(double));
        cudaMalloc((void**) &urecvDevice, u_recv.size() * sizeof(double));
        cudaMalloc((void**) &recvNeighboursDevice, recvNeighbours.size() * sizeof(Volume));

        cudaMemcpy(u0Device, u0.data(), u0.size() * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(u1Device, u1.data(), u1.size() * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(urecvDevice, u_recv.data(), u_recv.size() * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(recvNeighboursDevice, recvNeighbours.data(), recvNeighbours.size() * sizeof(Volume), cudaMemcpyHostToDevice);

        FillNextLayerKernel<<<nblocks, BLOCK_SIZE>>>(uDevice, u0Device, u1Device, urecvDevice, recvNeighboursDevice, total, recvNeighbours.size(), xmin, ymin, zmin, dy, dz, volume, hx, hy, hz, tau);
        cudaMemcpy(u.data(), uDevice, u.size() * sizeof(double), cudaMemcpyDeviceToHost);

        cudaFree(u0Device);
        cudaFree(u1Device);
        cudaFree(uDevice);
        cudaFree(urecvDevice);
        cudaFree(recvNeighboursDevice);
    }

    FillBoundaryValues(u, t, bt, hx, hy, hz, L, volume, N);
}

// заполнение аналитических значений
void FillAnalyticalValues(std::vector<double> &u, double t, double hx, double hy, double hz, VolumeSize L, Volume volume) {
    for (int i = volume.xmin; i <= volume.xmax; i++)
        for (int j = volume.ymin; j <= volume.ymax; j++)
            for (int k = volume.zmin; k <= volume.zmax; k++)
                u[Index(i, j, k, volume)] = AnalyticalSolve(i * hx, j * hy, k * hz, t, L);
}

// заполнение разности значений
void FillDifferenceValues(std::vector<double> &u, double t, double hx, double hy, double hz, VolumeSize L, Volume volume) {
    for (int i = volume.xmin; i <= volume.xmax; i++) {
        for (int j = volume.ymin; j <= volume.ymax; j++) {
            for (int k = volume.zmin; k <= volume.zmax; k++) {
                u[Index(i, j, k, volume)] = fabs(u[Index(i, j, k, volume)] - AnalyticalSolve(i * hx, j * hy, k * hz, t, L));
            }
        }
    }
}

// оценка погрешности на слое
double EvaluateError(const std::vector<double> &u, double t, double hx, double hy, double hz, VolumeSize L, Volume volume) {
    double localError = 0;

    if (!USE_GPU) {
        for (int i = volume.xmin; i <= volume.xmax; i++)
            for (int j = volume.ymin; j <= volume.ymax; j++)
                for (int k = volume.zmin; k <= volume.zmax; k++)
                    localError = std::max(localError, fabs(u[Index(i, j, k, volume)] - AnalyticalSolve(i * hx, j * hy, k * hz, t, L)));
    }
    else {
        double *uDevice, *errorDevice;
        cudaMalloc((void**) &uDevice, u.size() * sizeof(double));
        cudaMalloc((void**) &errorDevice, sizeof(double));
        cudaMemcpy(uDevice, u.data(), u.size() * sizeof(double), cudaMemcpyHostToDevice);

        EvaluateErrorKernel<<<1, BLOCK_SIZE>>>(uDevice, hx, hy, hz, t, L, volume, errorDevice);

        cudaMemcpy(&localError, errorDevice, 1 * sizeof(double), cudaMemcpyDeviceToHost);
        cudaFree(uDevice);
        cudaFree(errorDevice);
    }

    double error;
    MPI_Reduce(&localError, &error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return error;
}

// сохранение слоя
void SaveValues(const std::vector<double> u, double t, const std::vector<Volume> &volumes, int rank, int size, Volume volume, int layerSize, VolumeSize L, int N, const char *filename) {
    std::vector<double> u_all = SendRecvTotal(u, volumes, rank, size, volume, layerSize, N);

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
double CudaMPIHyperbolicEquationSolve(VolumeSize L, double T, int N, int K, BoundaryConditionTypes bt, SplitType split, int rank, int size, SolveParams params, bool isSilent) {
    double hx = L.x / N;
    double hy = L.y / N;
    double hz = L.z / N;
    double tau = T / K;

    int layerSize = (N + 1) * (N + 1) * (N + 1);

    GridSplitter splitter(N);
    std::vector<Volume> volumes = splitter.Split(split, size);
    Volume volume = volumes[rank];

    std::vector<int> processNeighbours;
    std::vector<Volume> sendNeighbours;
    std::vector<Volume> recvNeighbours;

    FillNeighbours(volumes, processNeighbours, sendNeighbours, recvNeighbours, rank, size, volume);

    std::vector<std::vector<double>> u(3, std::vector<double>(volume.size));

    FillInitialValues(u[0], u[1], bt, hx, hy, hz, tau, L, N, processNeighbours, sendNeighbours, recvNeighbours, volume);

    double error0 = EvaluateError(u[0], 0, hx, hy, hz, L, volume);
    double error1 = EvaluateError(u[1], tau, hx, hy, hz, L, volume);

    if (rank == 0 && !isSilent) {
        std::ofstream fout(params.outputPath ? params.outputPath : "output.txt", std::ios::app);
        fout << "Layer 0 max error: " << error0 << std::endl;
        fout << "Layer 1 max error: " << error1 << std::endl;
        fout.close();
    }

    for (int step = 2; step <= params.steps; step++) {
        FillNextLayer(u[(step + 1) % 3], u[(step + 2) % 3], u[step % 3], step * tau, hx, hy, hz, tau, L, N, bt, processNeighbours, sendNeighbours, recvNeighbours, volume);

        double error = EvaluateError(u[step % 3], step * tau, hx, hy, hz, L, volume);

        if (rank == 0 && !isSilent) {
            std::ofstream fout(params.outputPath ? params.outputPath : "output.txt", std::ios::app);
            fout << "Layer " << step << " max error: " << error << std::endl;
            fout.close();
        }
    }

    if (params.numericalPath) {
        SaveValues(u[params.steps % 3], params.steps * tau, volumes, rank, size, volume, layerSize, L, N, params.numericalPath);
    }

    if (params.differencePath) {
        FillDifferenceValues(u[params.steps % 3], params.steps * tau, hx, hy, hz, L, volume);
        SaveValues(u[params.steps % 3], params.steps * tau, volumes, rank, size, volume, layerSize, L, N, params.differencePath);
    }

    if (params.analyticalPath) {
        FillAnalyticalValues(u[0], params.steps * tau, hx, hy, hz, L, volume);
        SaveValues(u[0], params.steps * tau, volumes, rank, size, volume, layerSize, L, N, params.analyticalPath);
    }

    return EvaluateError(u[params.steps % 3], params.steps * tau, hx, hy, hz, L, volume);
}

// вывод параметров
void PrintParams(VolumeSize L, double T, int N, int K, BoundaryConditionTypes bt, SplitType split, int size, const char *outputPath) {
    std::ofstream fout(outputPath ? outputPath : "output.txt", std::ios::app);

    fout << "Processors: " << size << std::endl << std::endl;

    fout << "Lx: " << L.x << std::endl;
    fout << "Ly: " << L.y << std::endl;
    fout << "Lz: " << L.z << std::endl;
    fout << "T: " << T << std::endl << std::endl;

    fout << "N: " << N << std::endl;
    fout << "K: " << K << std::endl << std::endl;

    fout << "x_i = i*" << (L.x / N) << ", i = 0..." << N << std::endl;
    fout << "y_i = i*" << (L.y / N) << ", i = 0..." << N << std::endl;
    fout << "z_i = i*" << (L.z / N) << ", i = 0..." << N << std::endl;
    fout << "t_i = i*" << (T / K) << ", i = 0..." << K << std::endl << std::endl;

    fout << "Type of boundary condition along axis X: " << bt.x << std::endl;
    fout << "Type of boundary condition along axis Y: " << bt.y << std::endl;
    fout << "Type of boundary condition along axis Z: " << bt.z << std::endl << std::endl;

    fout << "Split strategy: " << split << std::endl << std::endl;
    fout.close();
}
