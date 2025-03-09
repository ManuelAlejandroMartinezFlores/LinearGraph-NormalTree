# Set the minimum version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name
project(NormalTree)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")


# Add executable (Your source file(s))
add_executable(Main src/main.cpp)
# add_executable(benchmark src/benchmark_test.cpp)

# If you're using an external library like nlohmann/json
# Find and include the nlohmann/json header if installed via Homebrew
find_package(nlohmann_json REQUIRED)
find_package(OpenMP REQUIRED)


# Link nlohmann_json to your project
target_link_libraries(Main PRIVATE nlohmann_json::nlohmann_json)
target_link_libraries(Main PRIVATE OpenMP::OpenMP_CXX)



