#include <iostream>
#include <fstream>
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

    HyperbolicEquationSolver solver(arguments.L, arguments.T, arguments.N, arguments.K, arguments.bt);

    if (arguments.debug) {
        solver.PrintParams(arguments.outputPath);
    }

    solver.Solve(arguments.steps, arguments.outputPath, arguments.numericalPath, arguments.analyticalPath, arguments.differencePath);
}