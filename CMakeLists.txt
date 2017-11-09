cmake_minimum_required (VERSION 2.6)
project (Genetics)

set(CMAKE_CXX_FLAGS "-O3 -W -Wall -pedantic -std=c++11")


add_executable (Genetics src/main.cpp)


enable_testing()
add_subdirectory(test/googletest)
include_directories(test/googletest/include test/googletest)


add_executable (Genetics_UnitTest test/main.cpp)

target_link_libraries(Genetics_UnitTest gtest gtest_main)
add_test(Genetics_UnitTest Genetics_UnitTest)
