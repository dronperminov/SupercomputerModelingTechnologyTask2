#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

class HyperbolicEquationSolver {
    double Lx, Ly, Lz; // параметры параллелепипеда
    double T; // время
    BoundaryConditionType btX, btY, btZ; // граничные условия
    int N; // размер пространственной сетки
    int K; // размер временной сетки
    double hx, hy, hz; // шаги пространственной сетки
    double tau; // шаг временной сетки

    double ***InitLayer(int ni, int nj, int nk) const; // выделение памяти под слой
    void FreeLayer(double ***layer, int ni, int nj, int nk) const; // освобождение памяти из-под слоя
    void SaveLayer(double ***layer, double t, const char *filename) const; // сохранение слоя

    void FillBoundaryValues(double ***u, double t, bool isInitial) const; // заполнение граничными значениями
    void FillInitialValues(double ***u0, double ***u1) const; // заполнение начальных условий

    double LaplaceOperator(double ***u, int i, int j, int k) const; // оператор Лапласа
    double EvaluateError(double ***u, double t) const; // оценка погрешности на слое
public:
    HyperbolicEquationSolver(double Lx, double Ly, double Lz, double T, int N, int K, BoundaryConditionType btX, BoundaryConditionType btY, BoundaryConditionType btZ);

    double AnalyticalSolve(double x, double y, double z, double t) const; // аналитическое решение
    double Phi(double x, double y, double z) const; // начальные условия

    void Solve(int maxSteps = 20, const char *jsonPath = NULL); // решение
    void PrintParams() const; // вывод параметров
};

HyperbolicEquationSolver::HyperbolicEquationSolver(double Lx, double Ly, double Lz, double T, int N, int K, BoundaryConditionType btX, BoundaryConditionType btY, BoundaryConditionType btZ) {
    this->Lx = Lx;
    this->Ly = Ly;
    this->Lz = Lz;
    this->T = T;

    this->N = N;
    this->K = K;

    this->hx = Lx / N;
    this->hy = Ly / N;
    this->hz = Lz / N;
    this->tau = T / K;

    this->btX = btX;
    this->btY = btY;
    this->btZ = btZ;
}

// аналитическое решение
double HyperbolicEquationSolver::AnalyticalSolve(double x, double y, double z, double t) const {
    double at = M_PI * sqrt(1 / (Lx * Lx) + 1 / (Ly * Ly) + 1 / (Lz * Lz));

    return sin(M_PI / Lx * x) * sin(2 * M_PI / Ly * y) * sin(3 * M_PI / Lz * z) * cos(at * t + 2*M_PI);
}

// начальные условия
double HyperbolicEquationSolver::Phi(double x, double y, double z) const {
    return AnalyticalSolve(x, y, z, 0);
}

// выделение памяти под слой
double*** HyperbolicEquationSolver::InitLayer(int ni, int nj, int nk) const {
    double ***layer = new double**[ni];

    for (int i = 0; i < ni; i++) {
        layer[i] = new double*[nj];

        for (int j = 0; j < nj; j++)
            layer[i][j] = new double[nk];
    }

    return layer;
}

// освобождение памяти из-под слоя
void HyperbolicEquationSolver::FreeLayer(double ***layer, int ni, int nj, int nk) const {
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++)
            delete[] layer[i][j];

        delete[] layer[i];
    }

    delete[] layer;
}

// сохранение слоя
void HyperbolicEquationSolver::SaveLayer(double ***layer, double t, const char *filename) const {
    std::ofstream f(filename);

    f << "{" << std::endl;
    f << "    \"Lx\": " << Lx << ", " << std::endl;
    f << "    \"Ly\": " << Ly << ", " << std::endl;
    f << "    \"Lz\": " << Lz << ", " << std::endl;
    f << "    \"N\": " << N << ", " << std::endl;
    f << "    \"t\": " << t << ", " << std::endl;
    f << "    \"points\": [" << std::endl;

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

                f << "        [" << (hx * i) << ", " << (hy * j) << ", " << (hz * k) << ", " << layer[i][j][k] << "]";
            }
        }
    }

    f << std::endl;
    f << "    ]" << std::endl;
    f << "}" << std::endl;

    f.close();
}

// оператор Лапласа
double HyperbolicEquationSolver::LaplaceOperator(double ***u, int i, int j, int k) const {
    double dx = (u[i - 1][j][k] - 2 * u[i][j][k] + u[i + 1][j][k]) / (hx * hx);
    double dy = (u[i][j - 1][k] - 2 * u[i][j][k] + u[i][j + 1][k]) / (hy * hy);
    double dz = (u[i][j][k - 1] - 2 * u[i][j][k] + u[i][j][k + 1]) / (hz * hz);

    return dx + dy + dz;
}

// заполнение граничными значениями
void HyperbolicEquationSolver::FillBoundaryValues(double ***u, double t, bool isInitial) const {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            if (btX == BoundaryConditionType::FirstKind) {
                u[0][i][j] = 0;
                u[N][i][j] = 0;
            }
            else if (btX == BoundaryConditionType::PeriodicAnalytical || (isInitial && btX == BoundaryConditionType::PeriodicNumerical)) {
                u[0][i][j] = AnalyticalSolve(0, i * hy, j * hz, t);
                u[N][i][j] = AnalyticalSolve(Lx, i * hy, j * hz, t);
            }
            else if (btX == BoundaryConditionType::PeriodicNumerical) {
                u[0][i][j] = (u[1][i][j] + u[N - 1][i][j]) / 2;
                u[N][i][j] = (u[1][i][j] + u[N - 1][i][j]) / 2;
            }

            if (btY == BoundaryConditionType::FirstKind) {
                u[i][0][j] = 0;
                u[i][N][j] = 0;
            }
            else if (btY == BoundaryConditionType::PeriodicAnalytical || (isInitial && btY == BoundaryConditionType::PeriodicNumerical)) {
                u[i][0][j] = AnalyticalSolve(i * hx, 0, j * hz, t);
                u[i][N][j] = AnalyticalSolve(i * hx, Ly, j * hz, t);
            }
            else if (btY == BoundaryConditionType::PeriodicNumerical) {
                u[i][0][j] = (u[i][1][j] + u[i][N - 1][j]) / 2;
                u[i][N][j] = (u[i][1][j] + u[i][N - 1][j]) / 2;
            }

            if (btZ == BoundaryConditionType::FirstKind) {
                u[i][j][0] = 0;
                u[i][j][N] = 0;
            }
            else if (btZ == BoundaryConditionType::PeriodicAnalytical || (isInitial && btZ == BoundaryConditionType::PeriodicNumerical)) {
                u[i][j][0] = AnalyticalSolve(i * hx, j * hy, 0, t);
                u[i][j][N] = AnalyticalSolve(i * hx, j * hy, Lz, t);
            }
            else if (btZ == BoundaryConditionType::PeriodicNumerical) {
                u[i][j][0] = (u[i][j][1] + u[i][j][N - 1]) / 2;
                u[i][j][N] = (u[i][j][1] + u[i][j][N - 1]) / 2;
            }
        }
    }
}

// заполнение начальных условий
void HyperbolicEquationSolver::FillInitialValues(double ***u0, double ***u1) const {
    FillBoundaryValues(u0, 0, true);
    FillBoundaryValues(u1, tau, true);
    
    #pragma omp parallel for collapse(3)
    for (int i = 1; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 1; k < N; k++)
                u0[i][j][k] = Phi(i * hx, j * hy, k * hz);

    #pragma omp parallel for collapse(3)
    for (int i = 1; i < N; i++)
        for (int j = 1; j < N; j++)
            for (int k = 1; k < N; k++)
                u1[i][j][k] = u0[i][j][k] + tau * tau / 2 * LaplaceOperator(u0, i, j, k);
}

// оценка погрешности на слое
double HyperbolicEquationSolver::EvaluateError(double ***u, double t) const {
    double error = 0;

    #pragma omp parallel for collapse(3) reduction(max: error)
    for (int i = 0; i <= N; i++)
        for (int j = 0; j <= N; j++)
            for (int k = 0; k <= N; k++)
                error = std::max(error, fabs(u[i][j][k] - AnalyticalSolve(i * hx, j * hy, k * hz, t)));

    return error;
}

// решение
void HyperbolicEquationSolver::Solve(int maxSteps, const char *jsonPath) {
    double ****u = new double***[maxSteps + 1];

    for (int i = 0; i <= maxSteps; i++)
        u[i] = InitLayer(N + 1, N + 1, N + 1);

    FillInitialValues(u[0], u[1]); // заполняем начальные условия

    std::cout << "Layer 0 max error: " << EvaluateError(u[0], 0) << std::endl;
    std::cout << "Layer 1 max error: " << EvaluateError(u[1], tau) << std::endl;

    for (int step = 2; step <= maxSteps; step++) {
        #pragma omp parallel for collapse(3)
        for (int i = 1; i < N; i++)
            for (int j = 1; j < N; j++)
                for (int k = 1; k < N; k++)
                    u[step][i][j][k] = 2 * u[step - 1][i][j][k] - u[step - 2][i][j][k] + tau * tau * LaplaceOperator(u[step - 1], i, j, k);

        FillBoundaryValues(u[step], step * tau, false);

        std::cout << "Layer " << step << " max error: " << EvaluateError(u[step], step * tau) << std::endl;
    }

    if (jsonPath) {
        SaveLayer(u[maxSteps], maxSteps * tau, jsonPath);
    }

    for (int i = 0; i <= maxSteps; i++)
        FreeLayer(u[i], N + 1, N + 1, N + 1);

    delete[] u;
}

// вывод параметров
void HyperbolicEquationSolver::PrintParams() const {
    std::cout << "Lx: " << Lx << std::endl;
    std::cout << "Ly: " << Ly << std::endl;
    std::cout << "Lz: " << Lz << std::endl;
    std::cout << "T: " << T << std::endl << std::endl;

    std::cout << "N: " << N << std::endl;
    std::cout << "K: " << K << std::endl << std::endl;

    std::cout << "x_i = i*" << hx << ", i = 0..." << N << std::endl;
    std::cout << "y_i = i*" << hy << ", i = 0..." << N << std::endl;
    std::cout << "z_i = i*" << hz << ", i = 0..." << N << std::endl;
    std::cout << "t_i = i*" << tau << ", i = 0..." << K << std::endl << std::endl;

    std::cout << "Type of boundary condition along axis X: " << btX << std::endl;
    std::cout << "Type of boundary condition along axis Y: " << btY << std::endl;
    std::cout << "Type of boundary condition along axis Z: " << btZ << std::endl;
}