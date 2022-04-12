//
// Created by Mackem Meya on 10.02.2022.
//

#include <iostream>
#include <cmath>

#include "WaveEquationSolver.h"
#include "Exceptions/InvalidIndexException.h"

double WaveEquationSolver::calculateImpulseSourceFunction(size_t step, size_t i, size_t j, size_t sourceX,
                                                          size_t sourceY) const {
    if (j != sourceX || i != sourceY) {
        return 0;
    }

    double expression{ 2.0 * M_PI * ((double) step * m_timeStep - 1.5) };

    return exp(-0.0625 * expression * expression) * sin(expression);
}

Grid WaveEquationSolver::calculatePhaseSpeed() const {
    Grid phaseSpeed{ m_gridHeight, m_gridWidth };

    size_t halfGridWidth{ size_t(0.5 * m_gridWidth) };
    size_t shift{ 0 };
    for (size_t i{ 0 }; i < m_gridHeight; ++i) {
        for (size_t j{ 0 }; j < halfGridWidth; ++j) {
            (*phaseSpeed)[shift + j] = 0.01;

        }

        for (size_t j{ size_t(halfGridWidth) }; j < m_gridWidth; ++j) {
            (*phaseSpeed)[shift + j] = 0.04;
        }
        shift += m_gridWidth;
    }

    return phaseSpeed;
}

Grid WaveEquationSolver::solve(size_t sourceX, size_t sourceY) {
    if (sourceX > m_gridWidth - 1) {
        throw InvalidIndexException{
                "Invalid source x: " + std::to_string(sourceX) + ". The value must be between 0 ... " +
                std::to_string(m_gridWidth - 1) };
    }

    if (sourceY > m_gridHeight - 1) {
        throw InvalidIndexException{
                "Invalid source y: " + std::to_string(sourceY) + ". The value must be between 0 ... " +
                std::to_string(m_gridHeight - 1) };
    }

    Grid phaseSpeed{ calculatePhaseSpeed() };
    Grid result{ m_gridHeight, m_gridWidth, 0.0 };
    Grid previousResult{ m_gridHeight, m_gridWidth, 0.0 };

    Grid *resultPtr{ &result };
    Grid *previousResultPtr{ &previousResult };

    double timeStepSquared{ m_timeStep * m_timeStep };
    double xGridStepModified{ 1 / (2 * m_xGridStep * m_xGridStep) };
    double yGridStepModified{ 1 / (2 * m_yGridStep * m_yGridStep) };

    double maxValue{ 0 };

    for (size_t step{ 0 }; step < m_numberOfSteps; ++step) {
        double impulseSourceValue{ calculateImpulseSourceFunction(step, sourceY, sourceX, sourceX, sourceY) };
        size_t shift{ m_gridWidth };
        for (size_t i{ 1 }; i < m_gridHeight - 1; ++i) {
            double *shiftedResult{ &(**resultPtr)[shift] };
            double *previousShiftedResult{ &(**previousResultPtr)[shift] };
            double *shiftedPhaseSpeed{ &(*phaseSpeed)[shift] };

            #pragma ivdep
            for (size_t j{ 1 }; j < m_gridWidth - 1; ++j) {
                double currentImpulseSourceValue{ 0 };
                if (i == sourceY && j == sourceX) {
                    currentImpulseSourceValue = impulseSourceValue;
                }

                double currentU{ previousShiftedResult[j] };
                double pij{ shiftedPhaseSpeed[j] };
                double pi1j{ shiftedPhaseSpeed[-m_gridWidth + j] };
                double pij1{ shiftedPhaseSpeed[j - 1] };
                double pi1j1{ shiftedPhaseSpeed[-m_gridWidth + j - 1] };

                double computedValue{ 2 * currentU - shiftedResult[j] + timeStepSquared *
                                                                        (currentImpulseSourceValue +
                                                                         ((previousShiftedResult[j + 1] - currentU) *
                                                                          (pi1j + pij) +
                                                                          (previousShiftedResult[j - 1] - currentU) *
                                                                          (pi1j1 + pij1)) *
                                                                         (xGridStepModified) +
                                                                         ((previousShiftedResult[m_gridWidth + j] -
                                                                           currentU) *
                                                                          (pij1 + pij) +
                                                                          (previousShiftedResult[-m_gridWidth + j] -
                                                                           currentU) *
                                                                          (pi1j1 + pi1j)) *
                                                                         (yGridStepModified)) };
                maxValue = std::max(std::abs(computedValue), maxValue);
                shiftedResult[j] = computedValue;
            }
            shift += m_gridWidth;
        }

        std::swap(resultPtr, previousResultPtr);
    }

    std::cout << "Max value = " << maxValue << std::endl;

    return result;
}
