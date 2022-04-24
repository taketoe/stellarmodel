cmake_minimum_required(VERSION 3.0.0)
project(stellar VERSION 0.1.0)

include(CTest)
enable_testing()

#link_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
include_directories("/opt/homebrew/Cellar/python@3.9/3.9.12/Frameworks/Python.framework/Versions/3.9/include/python3.9/")
include_directories("/opt/homebrew/lib/python3.9/site-packages/numpy/core/include/")
link_directories("/opt/python-3.10.4/lib")
link_libraries("python3.10")


add_executable(stellar src/main.cpp src/Stellar.cpp src/Stellar.h)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
