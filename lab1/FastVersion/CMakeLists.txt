cmake_minimum_required(VERSION 3.20)
project(WaveEquationSolver)

set(CMAKE_CXX_STANDARD 14)

add_executable(waveEquationSolver main.cpp Grid.h WaveEquationSolver.h FileWriter.h Exceptions/FileNotOpenedException.h Grid.cpp Exceptions/InvalidIndexException.h WaveEquationSolver.cpp Exceptions/BaseException.h Timer.cpp Timer.h Exceptions/MemoryNotAllocatedException.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")