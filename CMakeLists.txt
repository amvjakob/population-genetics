cmake_minimum_required (VERSION 2.6)
project (Genetics)

set(CMAKE_CXX_FLAGS "-O3 -W -Wall -pthread -pedantic -std=c++11")

set(SOURCE_FILES src/Allele.cpp src/Simulation.cpp src/SimulationsExecutor.cpp src/Random.cpp) 

add_executable (Genetics src/main.cpp ${SOURCE_FILES})

enable_testing()
add_subdirectory(test/gtest)
include_directories(test/gtest/include test/googletest)

add_executable (Genetics_UnitTest src/test_main.cpp src/random.cpp)
add_executable (testData src/data.cpp test/dataTest.cpp)

target_link_libraries(Genetics_UnitTest gtest gtest_main)
add_test(Genetics_UnitTest Genetics_UnitTest)

target_link_libraries(testData gtest gtest_main)
add_test(testData testData)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(tclap/include)

add_executable(program_using_tclap src/random_main.cpp src/random.cpp)

# We first check if Doxygen is present.
find_package(Doxygen)
if(DOXYGEN_FOUND)

    # We then check if the Doxygen.in configuration file exists
    IF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
        message(FATAL_ERROR "Doxygen configuration file is needed to build the documentation:
         - ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in not found.
         - Generate one in the project folder with doxygen -g Doxyfile.in")
    ENDIF(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)

    # We configure the doxyfile for make.
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
ELSE(DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
ENDIF(DOXYGEN_FOUND)
