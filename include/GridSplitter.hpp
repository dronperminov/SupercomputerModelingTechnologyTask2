#pragma once

#include <vector>
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

Volume MakeVolume(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax) {
    int dx = xmax - xmin + 1;
    int dy = ymax - ymin + 1;
    int dz = zmax - zmin + 1;
    int size = dx * dy * dz;

    return { xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, size };
}

class GridSplitter {
    int n; // размер сетки

    void SplitBlocks(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes);
    void SplitTapes(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes);
public:
    GridSplitter(int n);

    std::vector<Volume> Split(SplitType split, int size);
};

GridSplitter::GridSplitter(int n) {
    this->n = n;
}

void GridSplitter::SplitBlocks(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes) {
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
        SplitBlocks(xmin, x, ymin, ymax, zmin, zmax, size / 2, Y_AXIS, volumes);
        SplitBlocks(x + 1, xmax, ymin, ymax, zmin, zmax, size / 2, Y_AXIS, volumes);
    }
    else if (axis == Y_AXIS) {
        int y = (ymin + ymax) / 2;
        SplitBlocks(xmin, xmax, ymin, y, zmin, zmax, size / 2, Z_AXIS, volumes);
        SplitBlocks(xmin, xmax, y + 1, ymax, zmin, zmax, size / 2, Z_AXIS, volumes);
    }
    else {
        int z = (zmin + zmax) / 2;
        SplitBlocks(xmin, xmax, ymin, ymax, zmin, z, size / 2, X_AXIS, volumes);
        SplitBlocks(xmin, xmax, ymin, ymax, z + 1, zmax, size / 2, X_AXIS, volumes);
    }
}

void GridSplitter::SplitTapes(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int size, char axis, std::vector<Volume> &volumes) {
    int prev;

    if (axis == X_AXIS) {
        prev = xmin;
    }
    else if (axis == Y_AXIS) {
        prev = ymin;
    }
    else if (axis == Z_AXIS) {
        prev = zmin;
    }

    for (int i = 0; i < size; i++) {
        if (axis == X_AXIS) {
            int x = std::min(xmax, xmin + (xmax - xmin) * (i + 1) / size);
            volumes.push_back(MakeVolume(prev, x, ymin, ymax, zmin, zmax));
            prev = x + 1;
        }
        else if (axis == Y_AXIS) {
            int y = std::min(ymax, ymin + (ymax - ymin) * (i + 1) / size);
            volumes.push_back(MakeVolume(xmin, xmax, prev, y, zmin, zmax));
            prev = y + 1;
        }
        else if (axis == Z_AXIS) {
            int z = std::min(zmax, zmin + (zmax - zmin) * (i + 1) / size);
            volumes.push_back(MakeVolume(xmin, xmax, ymin, ymax, prev, z));
            prev = z + 1;
        }
    }
}

std::vector<Volume> GridSplitter::Split(SplitType split, int size) {
    std::vector<Volume> volumes;

    if (split == SplitType::Blocks) {
        SplitBlocks(0, n, 0, n, 0, n, size, X_AXIS, volumes);
    }
    else if (split == SplitType::Tapes) {
        SplitTapes(0, n, 0, n, 0, n, size, X_AXIS, volumes);
    }

    return volumes;
}