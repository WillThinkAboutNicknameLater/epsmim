cmake_minimum_required(VERSION 3.10)
project(CompilerVectorizedVersion)

set(CMAKE_CXX_COMPILER /glob/development-tools/versions/oneapi/2022.1.2/oneapi/compiler/2022.0.2/linux/bin/intel64/icc)
set(CMAKE_CXX_STANDARD 14)

add_executable(waveEquationSolver main.cpp Grid.h WaveEquationSolver.h FileWriter.h Exceptions/FileNotOpenedException.h Grid.cpp Exceptions/InvalidIndexException.h WaveEquationSolver.cpp Exceptions/BaseException.h Timer.cpp Timer.h Exceptions/MemoryNotAllocatedException.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")