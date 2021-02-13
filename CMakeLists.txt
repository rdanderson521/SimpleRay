cmake_minimum_required(VERSION 3.10)

project(SimpleRay VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

##### debug/dev options #####

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

##### libraries #####

find_package(X11 REQUIRED)
include_directories(${X11_INCLUDE_DIR})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(deps/CImg)

##### project code #####

set(sources src/main.cpp)
include_directories(src/)

add_executable(simpleray ${sources})

target_link_libraries(simpleray ${X11_LIBRARIES})
target_link_libraries(simpleray Threads::Threads)
