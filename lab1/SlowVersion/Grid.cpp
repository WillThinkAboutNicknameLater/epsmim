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
    for (size_t i{ 0 }; i < m_height; ++i) {
        for (size_t j{ 0 }; j < m_width; ++j) {
            std::cout << (*this)(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

void Grid::fill(double value) {
    for (size_t i{ 0 }; i < m_height; ++i) {
        for (size_t j{ 0 }; j < m_width; ++j) {
            (*this)(i, j) = value;
        }
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
        delete[] m_data;
        m_height = other.getHeight();
        m_width = other.getWidth();
        m_data = new double[m_height * m_width];
        if (m_data == nullptr) {
            throw MemoryNotAllocatedException{ "Failed to allocate memory for Grid" };
        }
    }

    for (size_t i{ 0 }; i < m_height; ++i) {
        for (size_t j{ 0 }; j < m_width; ++j) {
            (*this)(i, j) = other(i, j);
        }
    }

    return *this;
}