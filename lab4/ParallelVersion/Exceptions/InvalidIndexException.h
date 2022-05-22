//
// Created by Mackem Meya on 12.02.2022.
//

#ifndef INVALIDINDEXEXCEPTION_H
#define INVALIDINDEXEXCEPTION_H

#include "BaseException.h"

class InvalidIndexException : public BaseException {
public:
    explicit InvalidIndexException(const std::string &message) : BaseException{ message } { }
};

#endif //INVALIDINDEXEXCEPTION_H
