//
// Created by Mackem Meya on 12.02.2022.
//

#include <chrono>

#include "Timer.h"

double Timer::measure(const std::function<void()> &function) {
    auto start{ std::chrono::high_resolution_clock::now() };
    function();
    std::chrono::duration<double> elapsedTime{ std::chrono::high_resolution_clock::now() - start };
    return elapsedTime.count();
}
