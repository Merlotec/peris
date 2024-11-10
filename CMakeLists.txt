cmake_minimum_required(VERSION 3.18)
project(gspabm LANGUAGES CXX)

# Set the C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define the path to the locally installed SFML
set(SFML_DIR "C:/Users/ncbmk/dev/deps/SFML-2.6.1/lib/cmake/SFML")

# Find the SFML packages
find_package(SFML 2.6 COMPONENTS graphics window system REQUIRED)

# Add the executable or library
add_executable(gspabm
        src/main.cpp
        src/peris.cpp
        src/world.cpp "include/render.h" "include/solver.h")

add_definitions(-DPERIS_DISPLAY)


target_include_directories(gspabm PRIVATE include/)

# Link the SFML libraries
target_link_libraries(gspabm sfml-graphics sfml-window sfml-system)
