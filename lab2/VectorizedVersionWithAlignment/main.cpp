#include <iostream>

#include "WaveEquationSolver.h"
#include "Timer.h"
#include "FileWriter.h"
#include "Exceptions/InvalidParameterException.h"
#include "Exceptions/InvalidIndexException.h"
#include "Exceptions/FileNotOpenedException.h"

int main(int argc, char **argv) {
    if (argc != 7) {
        std::cerr << "Invalid arguments. Usage: \"" << argv[0]
                  << "\" gridWidth gridHeight numberOfSteps sourceX sourceY filePath" << std::endl;
        return 1;
    }

    size_t gridWidth{ (size_t) std::atol(argv[1]) };
    size_t gridHeight{ (size_t) std::atol(argv[2]) };
    size_t numberOfSteps{ (size_t) std::atol(argv[3]) };
    size_t sourceX{ (size_t) std::atol(argv[4]) };
    size_t sourceY{ (size_t) std::atol(argv[5]) };
    std::string filePath{ argv[6] };

    try {
        WaveEquationSolver waveEquationSolver{ gridWidth, gridHeight, numberOfSteps };
        Grid result{ };
        double elapsedTime{ Timer().measure([&result, &waveEquationSolver, &sourceX, &sourceY]() {
            result = std::move(waveEquationSolver.solve(sourceX, sourceY));
        }) };
        std::cout << "Elapsed time: " << elapsedTime << " s\n";

        FileWriter<double> fileWriter{ };
        fileWriter.write(filePath, *result, gridWidth * gridHeight);
    } catch (InvalidParameterException &e) {
        std::cerr << e.what() << std::endl;
    } catch (MemoryNotAllocatedException &e) {
        std::cerr << e.what() << std::endl;
    } catch (InvalidIndexException &e) {
        std::cerr << e.what() << std::endl;
    } catch (FileNotOpenedException &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}
