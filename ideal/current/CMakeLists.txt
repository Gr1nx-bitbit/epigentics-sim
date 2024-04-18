cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(peptides VERSION 1.0)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/aminoCodon.txt
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(GLOB_RECURSE SRC_FILES src/*.cpp)
add_executable(peptides main.cpp ${SRC_FILES})
target_include_directories(peptides PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)