#include <iostream>
#include <fstream>
#include <cmath>
#include "ArgumentParser.hpp"
#include "HyperbolicEquationSolver.hpp"

using namespace std;

int main(int argc, char **argv) {
    ArgumentParser parser;

    if (argc == 2 && (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))) {
        parser.Help();
        return 0;
    }

    Arguments arguments;

    try {
        arguments = parser.Parse(argc, argv);
    }
    catch (const char *error) {
        return -1;
    }

    HyperbolicEquationSolver solver(arguments.Lx, arguments.Ly, arguments.Lz, arguments.T, arguments.N, arguments.K, arguments.btX, arguments.btY, arguments.btZ);

    if (arguments.debug) {
        cout << "Readed parameters: " << endl;
        solver.PrintParams();
        cout << endl;
    }

    solver.Solve(arguments.steps, arguments.jsonPath);
}