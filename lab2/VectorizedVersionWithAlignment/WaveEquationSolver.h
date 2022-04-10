//
// Created by Mackem Meya on 10.02.2022.
//

#ifndef WAVEEQUATIONSOLVER_H
#define WAVEEQUATIONSOLVER_H

#include "Grid.h"

class WaveEquationSolver {
private:
    const size_t m_gridWidth;
    const size_t m_gridHeight;
    const size_t m_numberOfSteps;
    const double m_timeStep;
    const double m_xGridStep;
    const double m_yGridStep;

    double calculateImpulseSourceFunction(size_t step, size_t i, size_t j, size_t sourceX, size_t sourceY) const;

    Grid calculatePhaseSpeed() const;

public:
    WaveEquationSolver() = delete;

    WaveEquationSolver(size_t gridWidth, size_t gridHeight, size_t numberOfSteps) : m_gridWidth{ gridWidth },
                                                                                    m_gridHeight{ gridHeight },
                                                                                    m_numberOfSteps{ numberOfSteps },
                                                                                    m_timeStep{ gridWidth <= 1000 &&
                                                                                                gridHeight <= 1000
                                                                                                ? 0.01 : 0.001 },
                                                                                    m_xGridStep{ 4.0 /
                                                                                                 ((double) gridWidth -
                                                                                                  1.0) },
                                                                                    m_yGridStep{ 4.0 /
                                                                                                 ((double) gridHeight -
                                                                                                  1.0) } { }

    Grid solve(size_t sourceX, size_t sourceY);

    Grid solve2(size_t sourceX, size_t sourceY);

};

#endif //WAVEEQUATIONSOLVER_H
