//
// Created by Mackem Meya on 10.02.2022.
//

#include <iostream>
#include <cmath>
#include <immintrin.h>
#include <zmmintrin.h>

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
    __m256d timeStepSquaredV{ _mm256_set1_pd(timeStepSquared) };
    double xGridStepModified{ 1 / (2 * m_xGridStep * m_xGridStep) };
    __m256d xGridStepModifiedV{ _mm256_set1_pd(xGridStepModified) };
    double yGridStepModified{ 1 / (2 * m_yGridStep * m_yGridStep) };
    __m256d yGridStepModifiedV{ _mm256_set1_pd(yGridStepModified) };
    __m256d deucesV{ _mm256_set1_pd(2.0) };

    double maxValue{ 0 };
    __m256d maxValueV{ _mm256_setzero_pd() };
    __m256d maskForAbs{ _mm256_set1_pd(-0.) };

    for (size_t step{ 0 }; step < m_numberOfSteps; ++step) {
        double impulseSourceValue{ calculateImpulseSourceFunction(step, sourceY, sourceX, sourceX, sourceY) };
        size_t shift{ m_gridWidth };
        for (size_t i{ 1 }; i < m_gridHeight - 1; ++i) {

            double *shiftedResult{ &(**resultPtr)[shift] };
            double *previousShiftedResult{ &(**previousResultPtr)[shift] };
            double *shiftedPhaseSpeed{ &(*phaseSpeed)[shift] };

            __m256d *shiftedResultV{ (__m256d * )(shiftedResult + 1) };
            __m256d *left{ (__m256d * )(previousShiftedResult) };
            __m256d *center{ (__m256d * )(previousShiftedResult + 1) };
            __m256d *right{ (__m256d * )(previousShiftedResult + 2) };
            __m256d *leftP{ (__m256d * )(shiftedPhaseSpeed) };
            __m256d *down{ (__m256d * )(previousShiftedResult - m_gridWidth + 1) };
            __m256d *up{ (__m256d * )(previousShiftedResult + m_gridWidth + 1) };
            __m256d *downLeftP{ (__m256d * )(shiftedPhaseSpeed - m_gridWidth) };
            __m256d *centerP{ (__m256d * )(shiftedPhaseSpeed + 1) };
            __m256d *downCenterP{ (__m256d * )(shiftedPhaseSpeed - m_gridWidth + 1) };

            for (size_t j{ 1 }, vj{ 0 }; j < m_gridWidth - 4; j += 4, ++vj) {
                __m256d currentImpulseSourceValueV{ _mm256_setzero_pd() };
                if (i == sourceY && sourceX >= j && int(sourceX) - int(j) < 4) {
                    currentImpulseSourceValueV[sourceX - j] = impulseSourceValue;
                }
                __m256d currentCenter{ center[vj] };
                __m256d currentDownCenterP{ downCenterP[vj] };
                __m256d currentLeftP{ leftP[vj] };
                __m256d currentCenterP{ centerP[vj] };
                __m256d currentDownLeftP{ downLeftP[vj] };

                __m256d e1{ _mm256_fmsub_pd(deucesV, currentCenter, shiftedResultV[vj]) };
                __m256d e2{ (currentImpulseSourceValueV +
                             ((right[vj] - currentCenter) *
                              (currentDownCenterP +
                               currentCenterP) +
                              (left[vj] - currentCenter) *
                              (currentDownLeftP +
                               currentLeftP)) *
                             (xGridStepModifiedV) +
                             ((up[vj] -
                               currentCenter) *
                              (currentLeftP + currentCenterP) +
                              (down[vj] -
                               currentCenter) *
                              (currentDownLeftP +
                               currentDownCenterP)) *
                             (yGridStepModifiedV)) };

                __m256d computedValueV{ _mm256_fmadd_pd(timeStepSquaredV, e2, e1) };
                maxValueV = _mm256_max_pd(_mm256_andnot_pd(maskForAbs, computedValueV), maxValueV);
                shiftedResultV[vj] = computedValueV;
            }

            for (size_t j{ m_gridWidth - (m_gridWidth - 1) % 4 }; j < m_gridWidth - 1; ++j) {
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

    double *buffer{ (double *) (&maxValueV) };
    for (size_t i{ 0 }; i < 4; ++i) {
        maxValue = std::max(maxValue, buffer[i]);
    }
    std::cout << "Max value = " << maxValue << std::endl;

    return result;
}
