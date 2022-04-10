#include <iostream>

#include "WaveEquationSolver.h"
#include "Timer.h"
#include "FileWriter.h"
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

    WaveEquationSolver waveEquationSolver{ gridWidth, gridHeight, numberOfSteps };
    try {
        Grid result{ };
        double elapsedTime{ Timer().measure([&result, &waveEquationSolver, &sourceX, &sourceY]() {
            result = std::move(waveEquationSolver.solve(sourceX, sourceY));
        }) };
        std::cout << "Elapsed time: " << elapsedTime << " s\n";

        /*Grid result2{ };
        double elapsedTime2{ Timer().measure([&result2, &waveEquationSolver, &sourceX, &sourceY]() {
            result2 = std::move(waveEquationSolver.solve2(sourceX, sourceY));
        }) };
        std::cout << "Elapsed time: " << elapsedTime2 << " s\n";

        for (size_t i{ 0 }; i < result.getHeight(); ++i) {
            for (size_t j{ 0 }; j < result.getWidth(); ++j) {
                if (result(i, j) != result2(i, j)) {
                    std::cout << "False: i = " << i << ", j = " << j << std::endl;
                    std::cout << result(i, j) << " " << result2(i, j) << " " << result(i, j) - result2(i, j) << std::endl;
                    break;
                }
            }
        }*/

        FileWriter<double> fileWriter{ };
        fileWriter.write(filePath, *result, gridWidth * gridHeight);
    } catch (MemoryNotAllocatedException &e) {
        std::cerr << e.what() << std::endl;
    } catch (InvalidIndexException &e) {
        std::cerr << e.what() << std::endl;
    } catch (FileNotOpenedException &e) {
        std::cerr << e.what() << std::endl;
    }

    return 0;
}