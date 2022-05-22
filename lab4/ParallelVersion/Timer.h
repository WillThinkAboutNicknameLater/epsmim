//
// Created by Mackem Meya on 12.02.2022.
//

#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <functional>

class Timer {
public:
    double measure(const std::function<void()> &function);
};

#endif //TIMER_H
