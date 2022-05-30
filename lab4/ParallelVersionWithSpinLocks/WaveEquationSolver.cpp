//
// Created by Mackem Meya on 10.02.2022.
//

#include <iostream>
#include <cmath>
#include <immintrin.h>
#include <zmmintrin.h>
#include <omp.h>

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

void WaveEquationSolver::calculateRow(size_t rowIndex, Grid *resultPtr, Grid *previousResultPtr, Grid &phaseSpeed, size_t sourceX,
                                      size_t sourceY, double impulseSourceValue, double timeStepSquared, __m256d &timeStepSquaredV,
                                      double xGridStepModified, __m256d &xGridStepModifiedV,
                                      double yGridStepModified, __m256d &yGridStepModifiedV, __m256d &deucesV, __m256d &maskForAbs, double &maxValue,
                                      __m256d &maxValueV) const {
    size_t shift{ m_gridWidth * rowIndex };
    double *shiftedResult{ &(**resultPtr)[shift] };
    double *previousShiftedResult{ &(**previousResultPtr)[shift] };
    double *shiftedPhaseSpeed{ &(*phaseSpeed)[shift] };

    __m256d *shiftedResultV{ (__m256d *) (shiftedResult + 1) };

    __m256d *topU{ (__m256d *) (previousShiftedResult - m_gridWidth + 1) };
    __m256d *leftU{ (__m256d *) (previousShiftedResult) };
    __m256d *middleU{ (__m256d *) (previousShiftedResult + 1) };
    __m256d *rightU{ (__m256d *) (previousShiftedResult + 2) };
    __m256d *bottomU{ (__m256d *) (previousShiftedResult + m_gridWidth + 1) };

    __m256d *topLeftP{ (__m256d *) (shiftedPhaseSpeed - m_gridWidth) };
    __m256d *topP{ (__m256d *) (shiftedPhaseSpeed - m_gridWidth + 1) };
    __m256d *leftP{ (__m256d *) (shiftedPhaseSpeed) };
    __m256d *middleP{ (__m256d *) (shiftedPhaseSpeed + 1) };

    for (size_t j{ 1 }, vj{ 0 }; j < m_gridWidth - 4; j += 4, ++vj) {
        __m256d currentImpulseSourceValueV{ _mm256_setzero_pd() };
        if (rowIndex == sourceY && sourceX >= j && int(sourceX) - int(j) < 4) {
            currentImpulseSourceValueV[sourceX - j] = impulseSourceValue;
        }
        __m256d currentMiddleU{ middleU[vj] };
        __m256d currentTopP{ topP[vj] };
        __m256d currentLeftP{ leftP[vj] };
        __m256d currentMiddleP{ middleP[vj] };
        __m256d currentTopLeftP{ topLeftP[vj] };

        __m256d e1{ _mm256_fmsub_pd(deucesV, currentMiddleU, shiftedResultV[vj]) };
        __m256d e2{ (currentImpulseSourceValueV +
                     ((rightU[vj] - currentMiddleU) *
                      (currentTopP +
                       currentMiddleP) +
                      (leftU[vj] - currentMiddleU) *
                      (currentTopLeftP +
                       currentLeftP)) *
                     xGridStepModifiedV +
                     ((bottomU[vj] -
                       currentMiddleU) *
                      (currentLeftP + currentMiddleP) +
                      (topU[vj] -
                       currentMiddleU) *
                      (currentTopLeftP +
                       currentTopP)) *
                     yGridStepModifiedV) };

        __m256d computedValueV{ _mm256_fmadd_pd(timeStepSquaredV, e2, e1) };
        maxValueV = _mm256_max_pd(_mm256_andnot_pd(maskForAbs, computedValueV), maxValueV);
        shiftedResultV[vj] = computedValueV;
    }

    for (size_t j{ m_gridWidth - (m_gridWidth - 1) % 4 }; j < m_gridWidth - 1; ++j) {
        double currentImpulseSourceValue{ 0 };
        if (rowIndex == sourceY && j == sourceX) {
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

    double timeStepSquared{ m_timeStep * m_timeStep };
    __m256d timeStepSquaredV{ _mm256_set1_pd(timeStepSquared) };
    double xGridStepModified{ 1 / (2 * m_xGridStep * m_xGridStep) };
    __m256d xGridStepModifiedV{ _mm256_set1_pd(xGridStepModified) };
    double yGridStepModified{ 1 / (2 * m_yGridStep * m_yGridStep) };
    __m256d yGridStepModifiedV{ _mm256_set1_pd(yGridStepModified) };
    __m256d deucesV{ _mm256_set1_pd(2.0) };

    __m256d maskForAbs{ _mm256_set1_pd(-0.) };

    int numberOfThreads{ omp_get_max_threads() };

    size_t distance{ 16 };

    Grid *resultPtrPerThread[numberOfThreads * distance];
    Grid *previousResultPtrPerThread[numberOfThreads * distance];

    volatile size_t barriers[numberOfThreads * distance][2];

    size_t iterationsPerThread[numberOfThreads * distance];
    double maxValuePerThread[numberOfThreads * distance];
    __m256d maxValueVPerThread[numberOfThreads * distance];
    size_t intervalSize{ (m_gridHeight - 2) / numberOfThreads };
    for (size_t i{ 0 }; i < numberOfThreads; ++i) {
        resultPtrPerThread[i * distance] = &result;
        previousResultPtrPerThread[i * distance] = &previousResult;
        barriers[i * distance][0] = 0;
        barriers[i * distance][1] = 0;
        iterationsPerThread[i * distance] = intervalSize;
        maxValuePerThread[i * distance] = 0;
        maxValueVPerThread[i * distance] = _mm256_setzero_pd();
    }
    iterationsPerThread[(numberOfThreads - 1) * distance] += (m_gridHeight - 2) % numberOfThreads;

    #pragma omp parallel proc_bind(spread)
    for (size_t step{ 0 }; step < m_numberOfSteps; ++step) {
        int threadNumber{ omp_get_thread_num() };
        size_t shift{ threadNumber * distance };
        double impulseSourceValue{ calculateImpulseSourceFunction(step, sourceY, sourceX, sourceX, sourceY) };
        size_t startIndex{ 1 + threadNumber * intervalSize };
        size_t endIndex{ startIndex + iterationsPerThread[shift] - 1 };

        auto prevThreadBarrier{ barriers[(threadNumber - 1) * distance] };
        auto currentThreadBarrier{ barriers[shift] };
        auto nextThreadBarrier{ barriers[(threadNumber + 1) * distance] };

        if (threadNumber == 0) {
            for (size_t i{ startIndex }; i <= endIndex - 1; ++i) {
                calculateRow(i, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed, sourceX,
                             sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                             xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                             maxValuePerThread[shift], maxValueVPerThread[shift]);
            }
            while (currentThreadBarrier[1] > nextThreadBarrier[0]) {

            }
            calculateRow(endIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[1];
        } else if (threadNumber == numberOfThreads - 1) {
            while (currentThreadBarrier[0] > prevThreadBarrier[1]) {

            }
            calculateRow(startIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[0];

            for (size_t i{ startIndex + 1 }; i <= endIndex; ++i) {
                calculateRow(i, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed, sourceX,
                             sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                             xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                             maxValuePerThread[shift], maxValueVPerThread[shift]);
            }
        } else {
            bool flag{ false };
            size_t d{ 50 };
            size_t idk{ (endIndex - startIndex - 1) / d };
            for (size_t l{ 0 }; l < idk; ++l) {
                for (size_t i{ d * l + startIndex + 1 }; i < d * (l + 1) + startIndex + 1; ++i) {
                    calculateRow(i, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                                 sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                                 xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                                 maxValuePerThread[shift], maxValueVPerThread[shift]);
                }

                if (!flag && currentThreadBarrier[0] <= prevThreadBarrier[1]) {
                    calculateRow(startIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift],
                                 phaseSpeed, sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                                 xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                                 maxValuePerThread[shift], maxValueVPerThread[shift]);
                    ++currentThreadBarrier[0];
                    flag = true;
                }
            }

            for (size_t i{ d * idk + startIndex + 1 }; i <= endIndex - 1; ++i) {
                calculateRow(i, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed, sourceX,
                             sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                             xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                             maxValuePerThread[shift], maxValueVPerThread[shift]);
            }

            if (!flag) {
                while (currentThreadBarrier[0] > prevThreadBarrier[1]) {

                }
                calculateRow(startIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                             sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                             xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                             maxValuePerThread[shift], maxValueVPerThread[shift]);
                ++currentThreadBarrier[0];
            }

            while (currentThreadBarrier[1] > nextThreadBarrier[0]) {

            }
            calculateRow(endIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[1];
        }

        /*if (threadNumber == 0) {
            for (size_t i{ startIndex }; i <= endIndex - 1; ++i) {
                calculateRow(i, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed, sourceX,
                             sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                             xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                             maxValuePerThread[shift], maxValueVPerThread[shift]);
            }
            while (currentThreadBarrier[1] > nextThreadBarrier[0]) {

            }
            calculateRow(endIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[1];
        } else if (threadNumber == numberOfThreads - 1) {
            while (currentThreadBarrier[0] > prevThreadBarrier[1]) {

            }
            calculateRow(startIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[0];

            for (size_t i{ startIndex + 1 }; i <= endIndex; ++i) {
                calculateRow(i, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed, sourceX,
                             sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                             xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                             maxValuePerThread[shift], maxValueVPerThread[shift]);
            }
        } else {
            while (currentThreadBarrier[0] > prevThreadBarrier[1]) {

            }
            calculateRow(startIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[0];
            for (size_t i{ startIndex + 1 }; i <= endIndex - 1; ++i) {
                calculateRow(i, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed, sourceX,
                             sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                             xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                             maxValuePerThread[shift], maxValueVPerThread[shift]);
            }
            while (currentThreadBarrier[1] > nextThreadBarrier[0]) {

            }
            calculateRow(endIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[1];
        }*/

        /*for (size_t i{ startIndex + 1 }; i <= endIndex - 1; ++i) {
            calculateRow(i, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed, sourceX,
                         sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
        }
        if (threadNumber == 0) {
            calculateRow(startIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            while (currentThreadBarrier[1] > nextThreadBarrier[0]) {

            }
            calculateRow(endIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[1];
        } else if (threadNumber == numberOfThreads - 1) {
            calculateRow(endIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            while (currentThreadBarrier[0] > prevThreadBarrier[1]) {

            }
            calculateRow(startIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[0];
        } else {
            while (currentThreadBarrier[0] > prevThreadBarrier[1]) {

            }
            calculateRow(startIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[0];

            while (currentThreadBarrier[1] > nextThreadBarrier[0]) {

            }
            calculateRow(endIndex, resultPtrPerThread[shift], previousResultPtrPerThread[shift], phaseSpeed,
                         sourceX, sourceY, impulseSourceValue, timeStepSquared, timeStepSquaredV,
                         xGridStepModified, xGridStepModifiedV, yGridStepModified, yGridStepModifiedV, deucesV, maskForAbs,
                         maxValuePerThread[shift], maxValueVPerThread[shift]);
            ++currentThreadBarrier[1];
        }*/

        std::swap(resultPtrPerThread[shift], previousResultPtrPerThread[shift]);
    }

    double maxValue{ 0 };
    for (size_t i{ 0 }; i < numberOfThreads; ++i) {
        double *buffer{ (double *) (&maxValueVPerThread[i * distance]) };
        for (size_t j{ 0 }; j < 4; ++j) {
            maxValue = std::max(maxValue, buffer[j]);
        }
        maxValue = std::max(maxValue, maxValuePerThread[i * distance]);
    }

    std::cout << "Max value = " << maxValue << std::endl;

    return result;
}
