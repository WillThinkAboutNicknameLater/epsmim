//
// Created by Mackem Meya on 10.02.2022.
//

#ifndef WAVEEQUATIONSOLVER_H
#define WAVEEQUATIONSOLVER_H

#include <immintrin.h>
#include <zmmintrin.h>

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

    inline void calculateRow(size_t rowIndex, Grid *resultPtr, Grid *previousResultPtr, Grid &phaseSpeed, size_t sourceX,
                             size_t sourceY, double impulseSourceValue, double timeStepSquared, __m256d &timeStepSquaredV,
                             double xGridStepModified, __m256d &xGridStepModifiedV,
                             double yGridStepModified, __m256d &yGridStepModifiedV, __m256d &deucesV, __m256d &maskForAbs, double &maxValue,
                             __m256d &maxValueV) const;

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

};

#endif //WAVEEQUATIONSOLVER_H
