//
// Created by Mackem Meya on 12.02.2022.
//

#ifndef FILENOTOPENEDEXCEPTION_H
#define FILENOTOPENEDEXCEPTION_H

#include "BaseException.h"

class FileNotOpenedException : public BaseException {
public:
    explicit FileNotOpenedException(const std::string &message) : BaseException{ message } { }
};

#endif //FILENOTOPENEDEXCEPTION_H
