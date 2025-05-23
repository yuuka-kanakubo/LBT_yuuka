cmake_minimum_required (VERSION 3.16)
project (LBT CXX Fortran)


#install dependencies
install (DIRECTORY tables DESTINATION "${PROJECT_NAME}")
install (DIRECTORY hydroProfile DESTINATION "${PROJECT_NAME}")
install (FILES parameters.dat DESTINATION "${PROJECT_NAME}")
install (FILES parameters_g.dat parameters_ud.dat parameters_s.dat DESTINATION "${PROJECT_NAME}")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-Wall -Wextra -Wpedantic)
add_subdirectory (src)
