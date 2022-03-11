//
// Created by Mackem Meya on 12.02.2022.
//

#ifndef BASEEXCEPTION_H
#define BASEEXCEPTION_H

class BaseException : public std::exception {
protected:
    std::string m_message;

public:
    explicit BaseException(std::string message) : m_message{ std::move(message) } { }

    const char *what() const noexcept override { return m_message.c_str(); }
};

#endif //BASEEXCEPTION_H
