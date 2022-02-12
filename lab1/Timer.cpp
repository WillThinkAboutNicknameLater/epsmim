//
// Created by Mackem Meya on 12.02.2022.
//

#include "Timer.h"

std::chrono::duration<double> Timer::measure(const std::function<void()> &function) {
    auto start{ std::chrono::high_resolution_clock::now() };
    function();
    return std::chrono::high_resolution_clock::now() - start;
}
