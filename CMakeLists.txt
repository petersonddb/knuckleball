cmake_minimum_required(VERSION 3.28)
project(Knuckleball)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(nlohmann_json REQUIRED)

add_executable(knuckleball main.cpp)
target_link_libraries(knuckleball PRIVATE nlohmann_json::nlohmann_json)
