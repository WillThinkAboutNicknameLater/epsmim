cmake_minimum_required(VERSION 3.20)
project(lab1)

set(CMAKE_CXX_STANDARD 14)

add_executable(lab1 main.cpp Grid.h WaveEquationSolver.h FileWriter.h Exceptions/FileNotOpenedException.h Grid.cpp Exceptions/InvalidIndexException.h WaveEquationSolver.cpp Exceptions/BaseException.h Timer.cpp Timer.h)
