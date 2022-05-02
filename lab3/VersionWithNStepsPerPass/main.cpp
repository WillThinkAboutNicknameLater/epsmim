#include <iostream>

#include "WaveEquationSolver.h"
#include "Timer.h"
#include "FileWriter.h"
#include "Exceptions/InvalidParameterException.h"
#include "Exceptions/FileNotOpenedException.h"

int main(int argc, char **argv) {
    if (argc != 8) {
        std::cerr << "Invalid arguments. Usage: \"" << argv[0]
                  << "\" gridWidth gridHeight numberOfSteps sourceX sourceY numberOfStepsPerPass filePath" << std::endl;
        return 1;
    }

    size_t gridWidth{ (size_t) std::atol(argv[1]) };
    size_t gridHeight{ (size_t) std::atol(argv[2]) };
    size_t numberOfSteps{ (size_t) std::atol(argv[3]) };
    size_t sourceX{ (size_t) std::atol(argv[4]) };
    size_t sourceY{ (size_t) std::atol(argv[5]) };
    size_t numberOfStepsPerPass{ (size_t) std::atol(argv[6]) };
    std::string filePath{ argv[7] };

    WaveEquationSolver waveEquationSolver{ gridWidth, gridHeight, numberOfSteps };
    try {
        Grid result{ };
        double elapsedTime{ Timer().measure([&result, &waveEquationSolver, &sourceX, &sourceY, &numberOfStepsPerPass]() {
            result = std::move(waveEquationSolver.solve(sourceX, sourceY, numberOfStepsPerPass));
        }) };
        std::cout << "Elapsed time: " << elapsedTime << " s\n";

        FileWriter<double> fileWriter{ };
        fileWriter.write(filePath, *result, gridWidth * gridHeight);
    } catch (MemoryNotAllocatedException &e) {
        std::cerr << e.what() << std::endl;
    } catch (InvalidParameterException &e) {
        std::cerr << e.what() << std::endl;
    } catch (FileNotOpenedException &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}