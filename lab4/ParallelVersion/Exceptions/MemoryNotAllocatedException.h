//
// Created by Mackem Meya on 19.02.2022.
//

#ifndef MEMORYNOTALLOCATED_H
#define MEMORYNOTALLOCATED_H

#include "BaseException.h"

class MemoryNotAllocatedException : public BaseException {
public:
    explicit MemoryNotAllocatedException(const std::string &message) : BaseException{ message } { }
};

#endif //MEMORYNOTALLOCATED_H
