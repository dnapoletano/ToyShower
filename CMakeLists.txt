cmake_minimum_required(VERSION 3.0)
project(ToyShower++ LANGUAGES CXX)
set(CMAKE_MAKE_PROGRAM "Makefile")

message("------------ ToyShower++ --------------")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT RIVET_PREFIX)
  execute_process(COMMAND bash -c "rivet-config --prefix"
    OUTPUT_VARIABLE RIVET_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()


add_subdirectory(src)
message("------------   End Configure    --------------")