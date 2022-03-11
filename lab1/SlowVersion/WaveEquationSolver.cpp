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

    return exp(-(2.0 * M_PI * ((double) step * m_timeStep - 1.5)) * (2.0 * M_PI * ((double) step * m_timeStep - 1.5)) /
               (4.0 * 4.0)) *
           sin(2.0 * M_PI * ((double) step * m_timeStep - 1.5));
}

Grid WaveEquationSolver::calculatePhaseSpeed() const {
    Grid phaseSpeed{ m_gridHeight, m_gridWidth };

    for (size_t i{ 0 }; i < m_gridHeight; ++i) {
        for (size_t j{ 0 }; j < m_gridWidth / 2; ++j) {
            phaseSpeed(i, j) = 0.1 * 0.1;
        }

        for (size_t j{ m_gridWidth / 2 }; j < m_gridWidth; ++j) {
            phaseSpeed(i, j) = 0.2 * 0.2;
        }
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

    for (size_t step{ 0 }; step < m_numberOfSteps; ++step) {
        double impulseSourceValue{ calculateImpulseSourceFunction(step, sourceY, sourceX, sourceX, sourceY) };

        for (size_t i{ 1 }; i < m_gridHeight - 1; ++i) {
            for (size_t j{ 1 }; j < m_gridWidth - 1; ++j) {
                double currentImpulseSourceValue{ 0 };
                if (i == sourceY && j == sourceX) {
                    currentImpulseSourceValue = impulseSourceValue;
                }
                (*resultPtr)(i, j) = 2 * (*previousResultPtr)(i, j) - (*resultPtr)(i, j) + m_timeStep * m_timeStep *
                                                                                           (currentImpulseSourceValue +
                                                                                            (((*previousResultPtr)(i,
                                                                                                                   j +
                                                                                                                   1) -
                                                                                              (*previousResultPtr)(i,
                                                                                                                   j)) *
                                                                                             (phaseSpeed(i - 1, j) +
                                                                                              phaseSpeed(i, j)) +
                                                                                             ((*previousResultPtr)(i,
                                                                                                                   j -
                                                                                                                   1) -
                                                                                              (*previousResultPtr)(i,
                                                                                                                   j)) *
                                                                                             (phaseSpeed(i - 1, j - 1) +
                                                                                              phaseSpeed(i, j - 1))) /
                                                                                            (2 * m_xGridStep *
                                                                                             m_xGridStep) +
                                                                                            (((*previousResultPtr)(
                                                                                                    i + 1, j) -
                                                                                              (*previousResultPtr)(i,
                                                                                                                   j)) *
                                                                                             (phaseSpeed(i, j - 1) +
                                                                                              phaseSpeed(i, j)) +
                                                                                             ((*previousResultPtr)(
                                                                                                     i - 1, j) -
                                                                                              (*previousResultPtr)(i,
                                                                                                                   j)) *
                                                                                             (phaseSpeed(i - 1, j - 1) +
                                                                                              phaseSpeed(i - 1, j))) /
                                                                                            (2 * m_yGridStep *
                                                                                             m_yGridStep));
            }
        }

        std::swap(resultPtr, previousResultPtr);
    }

    return result;
}
