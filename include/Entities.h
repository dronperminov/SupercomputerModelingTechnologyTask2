#pragma once

// размер параллелепипеда
struct VolumeSize {
    double x; // длина параллелепипеда по x
    double y; // длина параллелепипеда по y
    double z; // длина параллелепипеда по z
};

// тип граничных условий
enum class BoundaryConditionType {
    FirstKind, // первого рода
    PeriodicAnalytical, // аналитические периодические
    PeriodicNumerical // численные периодические
};

// тип разбиения сетки
enum class SplitType {
    Blocks, // блочное разбиение
    Tapes, // ленточное разбиение
    Product // разбиение путём разложения на 3 сомножителя
};

// граничные условия
struct BoundaryConditionTypes {
    BoundaryConditionType x; // граничные условия по x
    BoundaryConditionType y; // граничные условия по y
    BoundaryConditionType z; // граничные условия по z
};

struct SolveParams {
    int steps; // количество шагов по времени
    char *outputPath; // путь для файла с выводом
    char *numericalPath; // путь для сохранения файла решения с точками
    char *analyticalPath; // путь для сохранения файла аналитического решения с точками
    char *differencePath; // путь для сохранения файла погрешности с точками
};

std::ostream& operator<<(std::ostream& os, const BoundaryConditionType type) {
    if (type == BoundaryConditionType::FirstKind)
        return os << "first-kind";

    if (type == BoundaryConditionType::PeriodicAnalytical)
        return os << "periodic-analytical";

    if (type == BoundaryConditionType::PeriodicNumerical)
        return os << "periodic-numerical";

    return os;
}

std::ostream& operator<<(std::ostream& os, const SplitType split) {
    if (split == SplitType::Blocks)
        return os << "blocks";

    if (split == SplitType::Tapes)
        return os << "tapes";

    if (split == SplitType::Product)
        return os << "product";

    return os;
}