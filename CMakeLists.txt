cmake_minimum_required (VERSION 2.6)
project (Genetics)

set(CMAKE_CXX_FLAGS "-O3 -W -Wall -pthread -pedantic -std=c++11")


set(SOURCE_FILES src/Allele.cpp src/Simulation.cpp src/SimulationsExecutor.cpp src/Random.cpp src/Data.cpp src/Globals.hpp) 

# Main executable
add_executable (Genetics src/main.cpp ${SOURCE_FILES})

# Test
enable_testing()
add_subdirectory(test/googletest)
include_directories(test/googletest/include test/googletest)
include_directories(src/tclap/include)

add_executable(testGenetics test/main.cpp ${SOURCE_FILES})

target_link_libraries(testGenetics gtest gtest_main)
add_test(testGenetics test)


###### Doxygen generation ######

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
