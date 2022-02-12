//
// Created by Mackem Meya on 08.02.2022.
//

#ifndef GRID_H
#define GRID_H

class Grid {
private:
    size_t m_height;
    size_t m_width;
    double *m_data;

public:
    Grid() = default;

    Grid(size_t gridHeight, size_t gridWidth, double initialValue = 0) : m_height{ gridHeight }, m_width{ gridWidth },
                                                                         m_data{ new double[gridHeight * gridWidth] } {
        fill(initialValue);
    }

    Grid(const Grid &other) : Grid() {
        *this = other;
    }

    ~Grid() {
        delete[] m_data;
    }

    size_t getHeight() const;

    size_t getWidth() const;

    void print() const;

    void fill(double value);

    double &operator()(size_t row, size_t column);

    const double &operator()(size_t row, size_t column) const;

    double *operator*();

    Grid& operator=(const Grid& other);

};

#endif //GRID_H
