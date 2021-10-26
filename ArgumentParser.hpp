#pragma once

#include <cstring>
#include <iostream>
#include "Entities.h"

struct Arguments {
    VolumeSize L; // размер параллелепипела
    double T; // итоговое время

    int N; // количество точек пространственной сетки
    int K; // количество точек временной сетки
    int steps; // количество шагов по времени

    BoundaryConditionTypes bt; // граничные условия

    bool debug; // режим отладки
    char *numericalPath; // путь для сохранения файла решения с точками
    char *analyticalPath; // путь для сохранения файла аналитического решения с точками
};

class ArgumentParser {
    bool IsInteger(const char *s) const;
    bool IsReal(const char *s) const;

    int ParseInteger(int argc, char **argv, int i, bool onlyPositive = false) const;
    double ParseReal(int argc, char **argv, int i, bool onlyPositive = false) const;

    BoundaryConditionType GetType(const char *arg) const; // получение типа граничного условия
public:
    void Help() const; // вывод сообщения помощи

    Arguments Parse(int argc, char **argv);
};

bool ArgumentParser::IsInteger(const char *s) const {
    for (int i = s[0] == '-' ? 1 : 0; s[i]; i++)
        if (s[i] < '0' || s[i] > '9')
            return false;

    return true;
}

bool ArgumentParser::IsReal(const char *s) const {
    int points = 0;
    int exponent = -1;
    int sign = -1;

    for (int i = s[0] == '-' ? 1 : 0; s[i]; i++) {
        if (s[i] == '.') {
            points++;

            if (points > 1 || exponent > -1)
                return false;
        }
        else if (s[i] == 'e') {
            if (exponent > -1)
                return false;

            exponent = i;
        }
        else if (s[i] == '-' || s[i] == '+') {
            if (exponent != i - 1)
                return false;

            sign = i;
        }
        else if (s[i] < '0' || s[i] > '9') {
            return false;
        }
    }

    if (exponent == -1)
        return true;

    if (sign == -1)
        return s[exponent + 1];

    return s[sign + 1];
}

int ArgumentParser::ParseInteger(int argc, char **argv, int i, bool onlyPositive) const {
    if (i >= argc - 1) {
        std::cout << argv[i] << " argument value missed" << std::endl;
        throw "argument value missed";
    }

    if (!IsInteger(argv[i + 1])) {
        std::cout << argv[i] << " argument value is not integer (" << argv[i + 1] << ")" << std::endl;
        throw "argument value is not integer";
    }

    int value = atoi(argv[i + 1]);

    if (onlyPositive && value <= 0) {
        std::cout << argv[i] << " argument value must be positive (" << argv[i + 1] << ")" << std::endl;
        throw "argument value must be positive";
    }

    return value;
}

double ArgumentParser::ParseReal(int argc, char **argv, int i, bool onlyPositive) const {
    if (i >= argc - 1) {
        std::cout << argv[i] << " argument value missed" << std::endl;
        throw "argument value missed";
    }

    if (!IsReal(argv[i + 1])) {
        std::cout << argv[i] << " argument value is not real (" << argv[i + 1] << ")" << std::endl;
        throw "argument value is not real";
    }

    double value = atof(argv[i + 1]);

    if (onlyPositive && value <= 0) {
        std::cout << argv[i] << " argument value must be positive (" << argv[i + 1] << ")" << std::endl;
        throw "argument value must be positive";
    }

    return value;
}

// получение типа граничного условия
BoundaryConditionType ArgumentParser::GetType(const char *arg) const {
    if (!strcmp(arg, "first") || !strcmp(arg, "first-kind") || !strcmp(arg, "f"))
        return BoundaryConditionType::FirstKind;

    if (!strcmp(arg, "periodic-analytical") || !strcmp(arg, "pa"))
        return BoundaryConditionType::PeriodicAnalytical;

    if (!strcmp(arg, "periodic-numerical") || !strcmp(arg, "pn"))
        return BoundaryConditionType::PeriodicNumerical;

    std::cout << "Unknown border condition type '" << arg << "'" << std::endl;
    throw "Unknown border condition type";
}

// вывод сообщения помощи
void ArgumentParser::Help() const {
    std::cout << "Usage: ./solve [-d] [-Lx Lx] [-Ly Ly] [-Lz Lz] [-T T] [-N N] [-K K] [-steps steps]" << std::endl;
    std::cout << "               [-btx border_type] [-bty border_type] [-btz border_type] [-oa path] [-on path]" << std::endl << std::endl;

    std::cout << "Arguments:" << std::endl;
    std::cout << "-d     - debug mode (default = non used)" << std::endl;
    std::cout << "-Lx    - the length of the parallelepiped along the X axis (default = 1)" << std::endl;
    std::cout << "-Ly    - the length of the parallelepiped along the Y axis (default = 1)" << std::endl;
    std::cout << "-Lz    - the length of the parallelepiped along the Z axis (default = 1)" << std::endl;
    std::cout << "-T     - the end time (default = 1)" << std::endl;
    std::cout << "-N     - number of points in dimension grid (default = 40)" << std::endl;
    std::cout << "-K     - number of points in time grid (default = 100)" << std::endl;
    std::cout << "-steps - number of steps for solving (default = 20)" << std::endl;
    std::cout << "-btx   - type of border condition along the X axis (default = first-kind)" << std::endl;
    std::cout << "-bty   - type of border condition along the Y axis (default = periodic-numerical)" << std::endl;
    std::cout << "-btz   - type of border condition along the Z axis (default = first-kind)" << std::endl;
    std::cout << "-on    - path to json file for saving numerical solve (default = non used)" << std::endl;
    std::cout << "-oa    - path to json file for saving analytical solve (default = non used)" << std::endl;
    std::cout << std::endl;
    std::cout << "Boundary condition types:" << std::endl;
    std::cout << "* first-kind (f)            - homogeneous boundary conditions of the first kind" << std::endl;
    std::cout << "* periodic-analytical (pa)  - analytic periodic boundary conditions" << std::endl;
    std::cout << "* periodic-numerical (pn)   - numerical periodic boundary conditions" << std::endl;
}

Arguments ArgumentParser::Parse(int argc, char **argv) {
    Arguments arguments;

    arguments.L.x = 1;
    arguments.L.y = 1;
    arguments.L.z = 1;
    arguments.T = 1;

    arguments.N = 40;
    arguments.K = 100;
    arguments.steps = 20;

    arguments.bt.x = BoundaryConditionType::FirstKind;
    arguments.bt.y = BoundaryConditionType::PeriodicNumerical;
    arguments.bt.z = BoundaryConditionType::FirstKind;

    arguments.debug = false;
    arguments.numericalPath = NULL;
    arguments.analyticalPath = NULL;

    for (int i = 1; i < argc; i += 2) {
        if (!strcmp(argv[i], "-Lx")) {
            arguments.L.x = ParseReal(argc, argv, i, true);
        }
        else if (!strcmp(argv[i], "-Ly")) {
            arguments.L.y = ParseReal(argc, argv, i, true);
        }
        else if (!strcmp(argv[i], "-Lz")) {
            arguments.L.z = ParseReal(argc, argv, i, true);
        }
        else if (!strcmp(argv[i], "-T")) {
            arguments.T = ParseReal(argc, argv, i);
        }
        else if (!strcmp(argv[i], "-N")) {
            arguments.N = ParseInteger(argc, argv, i, true);
        }
        else if (!strcmp(argv[i], "-K")) {
            arguments.K = ParseInteger(argc, argv, i, true);
        }
        else if (!strcmp(argv[i], "-steps")) {
            arguments.steps = ParseInteger(argc, argv, i, true);
        }
        else if (!strcmp(argv[i], "-btx")) {
            arguments.bt.x = GetType(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-bty")) {
            arguments.bt.y = GetType(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-btz")) {
            arguments.bt.z = GetType(argv[i + 1]);
        }
        else if (!strcmp(argv[i], "-on")) {
            if (i >= argc - 1)
                throw "-on argument value missed";

            arguments.numericalPath = argv[i + 1];
        }
        else if (!strcmp(argv[i], "-oa")) {
            if (i >= argc - 1)
                throw "-oa argument value missed";

            arguments.analyticalPath = argv[i + 1];
        }
        else if (!strcmp(argv[i], "-d") || !strcmp(argv[i], "--debug")) {
            arguments.debug = true;
            i--;
        }
        else {
            std::cout << "Unexpected argument '" << argv[i] << "'" << std::endl;
            throw "Unexpected argument";
        }
    }

    return arguments;
}