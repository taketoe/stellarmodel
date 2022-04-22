cmake_minimum_required(VERSION 3.0.0)
project(stellar VERSION 0.1.0)

include(CTest)
enable_testing()

#link_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")
include_directories("/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3")

add_executable(stellar src/main.cpp src/Stellar.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
