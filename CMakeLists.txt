cmake_minimum_required (VERSION 2.6)
project (LBT CXX Fortran)

add_subdirectory (src)

#install dependencies
install (DIRECTORY tables DESTINATION ${PROJECT_NAME})
install (DIRECTORY hydroProfile DESTINATION ${PROJECT_NAME})
install (FILES parameters.dat DESTINATION ${PROJECT_NAME})
#install (FILES parameters_g.dat parameters_ud.dat parameters_s.dat DESTINATION ${PROJECT_NAME})

