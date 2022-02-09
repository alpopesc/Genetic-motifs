cmake_minimum_required(VERSION 2.6)
project(Team8)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
SET(CMAKE_CXX_STRANDARD 11)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -pthread")
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -W -Wall -Wextra")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
option(test "Build tests." ON)
link_directories(${CMAKE_SOURCE_DIR}/lib)


set(SOURCE_FILES
        src/pwm.cpp
        src/motifs.cpp src/reader.cpp src/simulation.cpp
        src/extension.cpp
        )

add_executable(Team8 src/main.cpp ${SOURCE_FILES})
#SET(-DCMAKE_CXX_FLAGS=-D__GLIBCXX_USE_CXX11_ABI=0)



if (test)
    enable_testing()
    find_package(GTest)
    if (NOT GTEST_FOUND)
        SET(GTEST_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)
        SET(GTEST_BOTH_LIBRARIES libgtest.a libgtest_main.a)
    endif(NOT GTEST_FOUND)

    add_executable(testTeam8 tests/test_pwm.cpp tests/test_motifs.cpp tests/motifs_fasta_test.cpp
            tests/treader.cpp tests/test_extension.cpp ${SOURCE_FILES} tests/test_main.cpp)
    target_link_libraries(testTeam8 ${GTEST_BOTH_LIBRARIES} pthread)
    add_test(testTeam8 testTeam8)

endif(test)

find_package(Doxygen) 
if (DOXYGEN_FOUND)   
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)
