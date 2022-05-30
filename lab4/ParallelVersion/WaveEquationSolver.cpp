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

    __m256d maskForAbs{ _mm256_set1_pd(-0.) };

    int numberOfThreads{ omp_get_max_threads() };
    double maxValuePerThread[numberOfThreads];
    __m256d maxValueVPerThread[numberOfThreads];
    for (size_t i{ 0 }; i < numberOfThreads; ++i) {
        maxValuePerThread[i] = 0;
        maxValueVPerThread[i] = _mm256_setzero_pd();
    }

    for (size_t step{ 0 }; step < m_numberOfSteps; ++step) {
        double impulseSourceValue{ calculateImpulseSourceFunction(step, sourceY, sourceX, sourceX, sourceY) };

        #pragma omp parallel for proc_bind(spread)
        for (size_t i{ 1 }; i < m_gridHeight - 1; ++i) {
            int threadNumber{ omp_get_thread_num() };
            size_t shift{ m_gridWidth * i };
            double *shiftedResult{ &(**resultPtr)[shift] };
            double *previousShiftedResult{ &(**previousResultPtr)[shift] };
            double *shiftedPhaseSpeed{ &(*phaseSpeed)[shift] };

            __m256d *shiftedResultV{ (__m256d *) (shiftedResult + 1) };

            __m256d *topU{ (__m256d *) (previousShiftedResult - m_gridWidth) };
            __m256d *middleU{ (__m256d *) (previousShiftedResult) };
            __m256d *bottomU{ (__m256d *) (previousShiftedResult + m_gridWidth) };

            __m256d *topP{ (__m256d *) (shiftedPhaseSpeed - m_gridWidth) };
            __m256d *middleP{ (__m256d *) (shiftedPhaseSpeed) };

            __m256d topFirstU{ topU[0] };
            __m256d firstU{ middleU[0] };
            __m256d bottomFirstU{ bottomU[0] };

            __m256d topFirstP{ topP[0] };
            __m256d firstP{ middleP[0] };

            for (size_t j{ 1 }, vj{ 0 }; j < m_gridWidth - 4; j += 4, ++vj) {
                __m256d currentImpulseSourceValueV{ _mm256_setzero_pd() };
                if (i == sourceY && sourceX >= j && int(sourceX) - int(j) < 4) {
                    currentImpulseSourceValueV = _mm256_maskz_mov_pd(1 << (sourceX - j), _mm256_set1_pd(impulseSourceValue));
                }

                __m256d topSecondU{ topU[vj + 1] };
                __m256d secondU{ middleU[vj + 1] };
                __m256d bottomSecondU{ bottomU[vj + 1] };

                __m256d topSecondP{ topP[vj + 1] };
                __m256d secondP{ middleP[vj + 1] };

                __m256d topCenterU{ _mm256_permute4x64_pd(_mm256_blend_pd(topFirstU, topSecondU, 0b0001), 0b00111001) };
                __m256d centerU{ _mm256_permute4x64_pd(_mm256_blend_pd(firstU, secondU, 0b0001), 0b00111001) };
                __m256d bottomCenterU{ _mm256_permute4x64_pd(_mm256_blend_pd(bottomFirstU, bottomSecondU, 0b0001), 0b00111001) };

                __m256d rightU{ _mm256_permute4x64_pd(_mm256_blend_pd(firstU, secondU, 0b0011), 0b01001110) };

                __m256d topCenterP{ _mm256_permute4x64_pd(_mm256_blend_pd(topFirstP, topSecondP, 0b0001), 0b00111001) };
                __m256d centerP{ _mm256_permute4x64_pd(_mm256_blend_pd(firstP, secondP, 0b0001), 0b00111001) };

                __m256d e1{ _mm256_fmsub_pd(deucesV, centerU, shiftedResultV[vj]) };
                __m256d e2{ (currentImpulseSourceValueV +
                             ((rightU - centerU) *
                              (topCenterP + centerP) +
                              (firstU - centerU) *
                              (topFirstP + firstP)) *
                             xGridStepModifiedV +
                             ((bottomCenterU - centerU) *
                              (firstP + centerP) +
                              (topCenterU - centerU) *
                              (topFirstP + topCenterP)) *
                             yGridStepModifiedV) };

                __m256d computedValueV{ _mm256_fmadd_pd(timeStepSquaredV, e2, e1) };
                maxValueVPerThread[threadNumber] = _mm256_max_pd(_mm256_andnot_pd(maskForAbs, computedValueV), maxValueVPerThread[threadNumber]);
                shiftedResultV[vj] = computedValueV;

                topFirstU = topSecondU;
                firstU = secondU;
                bottomFirstU = bottomSecondU;

                topFirstP = topSecondP;
                firstP = secondP;
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
                                                                         xGridStepModified +
                                                                         ((previousShiftedResult[m_gridWidth + j] -
                                                                           currentU) *
                                                                          (pij1 + pij) +
                                                                          (previousShiftedResult[-m_gridWidth + j] -
                                                                           currentU) *
                                                                          (pi1j1 + pi1j)) *
                                                                         yGridStepModified) };

                maxValuePerThread[threadNumber] = std::max(std::abs(computedValue), maxValuePerThread[threadNumber]);
                shiftedResult[j] = computedValue;
            }
            shift += m_gridWidth;
        }

        std::swap(resultPtr, previousResultPtr);
    }

    double maxValue{ 0 };
    for (size_t i{ 0 }; i < numberOfThreads; ++i) {
        double *buffer{ (double *) (&maxValueVPerThread[i]) };
        for (size_t j{ 0 }; j < 4; ++j) {
            maxValue = std::max(maxValue, buffer[j]);
        }
        maxValue = std::max(maxValue, maxValuePerThread[i]);
    }
    
    std::cout << "Max value = " << maxValue << std::endl;

    return result;
}
