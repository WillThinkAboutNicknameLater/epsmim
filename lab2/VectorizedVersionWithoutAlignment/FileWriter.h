//
// Created by Mackem Meya on 10.02.2022.
//

#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <fstream>

#include "Exceptions/FileNotOpenedException.h"

template<typename T>
class FileWriter {
public:
    void write(const std::string &filePath, T *data, size_t size);
};

template<typename T>
void FileWriter<T>::write(const std::string &filePath, T *data, size_t size) {
    std::ofstream file{ filePath, std::ios::binary };

    if (!file.is_open()) {
        throw FileNotOpenedException{ "Failed to open file: " + filePath };
    }

    file.write((char *) data, size * sizeof(T));
    file.close();
}

#endif //FILEWRITER_H
