//
// Created by Mackem Meya on 11.04.2022.
//

#ifndef INVALIDPARAMETEREXCEPTION_H
#define INVALIDPARAMETEREXCEPTION_H

#include "BaseException.h"

class InvalidParameterException : public BaseException {
public:
    explicit InvalidParameterException(const std::string &message) : BaseException{ message } { }
};

#endif //INVALIDPARAMETEREXCEPTION_H
