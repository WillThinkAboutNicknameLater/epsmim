//
// Created by Mackem Meya on 08.02.2022.
//

#include <iostream>

#include "Grid.h"
#include "Exceptions/InvalidIndexException.h"

size_t Grid::getHeight() const {
    return m_height;
}

size_t Grid::getWidth() const {
    return m_width;
}

void Grid::print() const {
    size_t shift{ 0 };
    for (size_t i{ 0 }; i < m_height; ++i) {
        for (size_t j{ 0 }; j < m_width; ++j) {
            std::cout << m_data[shift + j] << " ";
        }
        shift += m_width;
        std::cout << std::endl;
    }
}

void Grid::fill(double value) {
    size_t shift{ 0 };
    for (size_t i{ 0 }; i < m_height; ++i) {
        for (size_t j{ 0 }; j < m_width; ++j) {
            m_data[shift + j] = value;
        }
        shift += m_width;
    }
}

double &Grid::operator()(size_t row, size_t column) {
    if (row >= m_height) {
        throw InvalidIndexException{
                "Invalid row number: " + std::to_string(row) + ". The value must be between 0 ... " +
                std::to_string(m_height - 1) };
    }

    if (column >= m_width) {
        throw InvalidIndexException{
                "Invalid column number: " + std::to_string(column) + ". The value must be between 0 ... " +
                std::to_string(m_width - 1) };
    }

    return m_data[row * m_width + column];
}

const double &Grid::operator()(size_t row, size_t column) const {
    if (row >= m_height) {
        throw InvalidIndexException{
                "Invalid row number: " + std::to_string(row) + ". The value must be between 0 ... " +
                std::to_string(m_height - 1) };
    }

    if (column >= m_width) {
        throw InvalidIndexException{
                "Invalid column number: " + std::to_string(column) + ". The value must be between 0 ... " +
                std::to_string(m_width - 1) };
    }

    return m_data[row * m_width + column];
}

double *Grid::operator*() {
    return m_data;
}

const double *Grid::operator*() const {
    return m_data;
}

Grid &Grid::operator=(const Grid &other) {
    if (this == &other) {
        return *this;
    }

    if (m_height != other.getHeight() && m_width != other.getWidth()) {
        free(m_data);
        m_height = other.getHeight();
        m_width = other.getWidth();
        m_data = (double *) aligned_alloc(32, m_height * m_width * 8);
        if (m_data == nullptr) {
            throw MemoryNotAllocatedException{ "Failed to allocate memory for Grid" };
        }
    }

    size_t shift{ 0 };
    for (size_t i{ 0 }; i < m_height; ++i) {
        for (size_t j{ 0 }; j < m_width; ++j) {
            m_data[shift + j] = (*other)[shift + j];
        }
        shift += m_width;
    }

    return *this;
}

Grid &Grid::operator=(Grid &&other) noexcept {
    if (this == &other) {
        return *this;
    }

    free(m_data);

    m_height = other.m_height;
    m_width = other.m_width;
    m_data = other.m_data;
    other.m_height = 0;
    other.m_width = 0;
    other.m_data = nullptr;

    return *this;
}
