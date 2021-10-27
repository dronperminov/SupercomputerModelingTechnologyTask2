#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include "Entities.h"

class HyperbolicEquationSolver {
    VolumeSize L; // параметры параллелепипеда
    double T; // время
    BoundaryConditionTypes bt; // граничные условия
    int N; // размер пространственной сетки
    int K; // размер временной сетки
    double hx, hy, hz; // шаги пространственной сетки
    double tau; // шаг временной сетки

    int layerSize; // размер слоя

    int Index(int i, int j, int k) const; // одномерная индексация
    void SaveLayer(double *layer, double t, const char *filename) const; // сохранение слоя

    void FillBoundaryValues(double *u, double t, bool isInitial) const; // заполнение граничными значениями
    void FillInitialValues(double *u0, double *u1) const; // заполнение начальных условий
    void FillAnalyticalValues(double *u, double t) const; // заполнение аналитических значений

    double LaplaceOperator(double *u, int i, int j, int k) const; // оператор Лапласа
    double EvaluateError(double *u, double t) const; // оценка погрешности на слое
public:
    HyperbolicEquationSolver(VolumeSize L, double T, int N, int K, BoundaryConditionTypes bt);

    double AnalyticalSolve(double x, double y, double z, double t) const; // аналитическое решение
    double Phi(double x, double y, double z) const; // начальные условия

    void Solve(int maxSteps = 20, const char *outputPath = NULL, const char *numericalPath = NULL, const char *analyticalPath = NULL); // решение
    void PrintParams(const char *outputPath) const; // вывод параметров
};

HyperbolicEquationSolver::HyperbolicEquationSolver(VolumeSize L, double T, int N, int K, BoundaryConditionTypes bt) {
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
}

// аналитическое решение
double HyperbolicEquationSolver::AnalyticalSolve(double x, double y, double z, double t) const {
    double at = M_PI * sqrt(1 / (L.x * L.x) + 1 / (L.y * L.y) + 1 / (L.z * L.z));

    return sin(M_PI / L.x * x) * sin(2 * M_PI / L.y * y) * sin(3 * M_PI / L.z * z) * cos(at * t + 2*M_PI);
}

// начальные условия
double HyperbolicEquationSolver::Phi(double x, double y, double z) const {
    return AnalyticalSolve(x, y, z, 0);
}

// одномерная индексация
int HyperbolicEquationSolver::Index(int i, int j, int k) const {
    return (i * (N + 1) + j) * (N + 1) + k;
}

// сохранение слоя
void HyperbolicEquationSolver::SaveLayer(double *layer, double t, const char *filename) const {
    std::ofstream f(filename);

    f << "{" << std::endl;
    f << "    \"Lx\": " << L.x << ", " << std::endl;
    f << "    \"Ly\": " << L.y << ", " << std::endl;
    f << "    \"Lz\": " << L.z << ", " << std::endl;
    f << "    \"N\": " << N << ", " << std::endl;
    f << "    \"t\": " << t << ", " << std::endl;
    f << "    \"u\": [" << std::endl;

    bool wasPrinted = false;

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            for (int k = 0; k <= N; k++) {
                if (wasPrinted) {
                    f << ", " << std::endl;
                }
                else {
                    wasPrinted = true;
                }

                f << "    " << layer[Index(i, j, k)];
            }
        }
    }

    f << std::endl;
    f << "    ]" << std::endl;
    f << "}" << std::endl;

    f.close();
}

// оператор Лапласа
double HyperbolicEquationSolver::LaplaceOperator(double *u, int i, int j, int k) const {
    double dx = (u[Index(i - 1, j, k)] - 2 * u[Index(i, j, k)] + u[Index(i + 1, j, k)]) / (hx * hx);
    double dy = (u[Index(i, j - 1, k)] - 2 * u[Index(i, j, k)] + u[Index(i, j + 1, k)]) / (hy * hy);
    double dz = (u[Index(i, j, k - 1)] - 2 * u[Index(i, j, k)] + u[Index(i, j, k + 1)]) / (hz * hz);

    return dx + dy + dz;
}

// заполнение граничными значениями
void HyperbolicEquationSolver::FillBoundaryValues(double *u, double t, bool isInitial) const {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            if (bt.x == BoundaryConditionType::FirstKind) {
                u[Index(0, i, j)] = 0;
                u[Index(N, i, j)] = 0;
            }
            else if (bt.x == BoundaryConditionType::PeriodicAnalytical || (isInitial && bt.x == BoundaryConditionType::PeriodicNumerical)) {
                u[Index(0, i, j)] = AnalyticalSolve(0, i * hy, j * hz, t);
                u[Index(N, i, j)] = AnalyticalSolve(L.x, i * hy, j * hz, t);
            }
            else if (bt.x == BoundaryConditionType::PeriodicNumerical) {
                u[Index(0, i, j)] = (u[Index(1, i, j)] + u[Index(N - 1, i, j)]) / 2;
                u[Index(N, i, j)] = (u[Index(1, i, j)] + u[Index(N - 1, i, j)]) / 2;
            }

            if (bt.y == BoundaryConditionType::FirstKind) {
                u[Index(i, 0, j)] = 0;
                u[Index(i, N, j)] = 0;
            }
            else if (bt.y == BoundaryConditionType::PeriodicAnalytical || (isInitial && bt.y == BoundaryConditionType::PeriodicNumerical)) {
                u[Index(i, 0, j)] = AnalyticalSolve(i * hx, 0, j * hz, t);
                u[Index(i, N, j)] = AnalyticalSolve(i * hx, L.y, j * hz, t);
            }
            else if (bt.y == BoundaryConditionType::PeriodicNumerical) {
                u[Index(i, 0, j)] = (u[Index(i, 1, j)] + u[Index(i, N - 1, j)]) / 2;
                u[Index(i, N, j)] = (u[Index(i, 1, j)] + u[Index(i, N - 1, j)]) / 2;
            }

            if (bt.z == BoundaryConditionType::FirstKind) {
                u[Index(i, j, 0)] = 0;
                u[Index(i, j, N)] = 0;
            }
            else if (bt.z == BoundaryConditionType::PeriodicAnalytical || (isInitial && bt.z == BoundaryConditionType::PeriodicNumerical)) {
                u[Index(i, j, 0)] = AnalyticalSolve(i * hx, j * hy, 0, t);
                u[Index(i, j, N)] = AnalyticalSolve(i * hx, j * hy, L.z, t);
            }
            else if (bt.z == BoundaryConditionType::PeriodicNumerical) {
                u[Index(i, j, 0)] = (u[Index(i, j, 1)] + u[Index(i, j, N - 1)]) / 2;
                u[Index(i, j, N)] = (u[Index(i, j, 1)] + u[Index(i, j, N - 1)]) / 2;
            }
        }
    }
}

// заполнение начальных условий
void HyperbolicEquationSolver::FillInitialValues(double *u0, double *u1) const {
    FillBoundaryValues(u0, 0, true);
    FillBoundaryValues(u1, tau, true);
    
    #pragma omp parallel for collapse(3)
    for (int i = 1; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 1; k < N; k++)
                u0[Index(i, j, k)] = Phi(i * hx, j * hy, k * hz);

    #pragma omp parallel for collapse(3)
    for (int i = 1; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 1; k < N; k++)
                u1[Index(i, j, k)] = u0[Index(i, j, k)] + tau * tau / 2 * LaplaceOperator(u0, i, j, k);
}

// заполнение аналитических значений
void HyperbolicEquationSolver::FillAnalyticalValues(double *u, double t) const {
    #pragma omp parallel for collapse(3)
    for (int i = 0; i <= N; i++)
        for (int j = 0; j <= N; j++)
            for (int k = 0; k <= N; k++)
                u[Index(i, j, k)] = AnalyticalSolve(i * hx, j * hy, k * hz, t);
}

// оценка погрешности на слое
double HyperbolicEquationSolver::EvaluateError(double *u, double t) const {
    double error = 0;

    #pragma omp parallel for collapse(3) reduction(max: error)
    for (int i = 0; i <= N; i++)
        for (int j = 0; j <= N; j++)
            for (int k = 0; k <= N; k++)
                error = std::max(error, fabs(u[Index(i, j, k)] - AnalyticalSolve(i * hx, j * hy, k * hz, t)));

    return error;
}

// решение
void HyperbolicEquationSolver::Solve(int maxSteps, const char *outputPath, const char *numericalPath, const char *analyticalPath) {
    double **u = new double*[3];
    u[0] = new double[layerSize];
    u[1] = new double[layerSize];
    u[2] = new double[layerSize];

    FillInitialValues(u[0], u[1]); // заполняем начальные условия

    std::ofstream fout(outputPath ? outputPath : "output.txt", std::ios::app);
    fout << "Layer 0 max error: " << EvaluateError(u[0], 0) << std::endl;
    fout << "Layer 1 max error: " << EvaluateError(u[1], tau) << std::endl;

    for (int step = 2; step <= maxSteps; step++) {
        #pragma omp parallel for collapse(3)
        for (int i = 1; i < N; i++)
            for (int j = 1; j < N; j++)
                for (int k = 1; k < N; k++)
                    u[step % 3][Index(i, j, k)] = 2 * u[(step + 2) % 3][Index(i, j, k)] - u[(step + 1) % 3][Index(i, j, k)] + tau * tau * LaplaceOperator(u[(step + 2) % 3], i, j, k);

        FillBoundaryValues(u[step % 3], step * tau, false);

        fout << "Layer " << step << " max error: " << EvaluateError(u[step % 3], step * tau) << std::endl;
    }

    fout.close();

    if (numericalPath) {
        SaveLayer(u[maxSteps % 3], maxSteps * tau, numericalPath);
    }

    if (analyticalPath) {
        FillAnalyticalValues(u[0], maxSteps * tau);
        SaveLayer(u[0], maxSteps * tau, analyticalPath);
    }

    delete[] u[0];
    delete[] u[1];
    delete[] u[2];
    delete[] u;
}

// вывод параметров
void HyperbolicEquationSolver::PrintParams(const char *outputPath) const {
    std::ofstream fout(outputPath ? outputPath : "output.txt", std::ios::app);

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
    fout << "Type of boundary condition along axis Z: " << bt.z << std::endl << std::endl;;

    fout.close();
}