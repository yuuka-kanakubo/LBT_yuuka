add_library(
		LBT_LIB STATIC
		LBTcl.cpp
	   )

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${CMAKE_HOME_DIRECTORY}
   )
add_executable(${PROJECT_NAME}
  main.cpp
  JetOutput-0.7.f
  hydro-CCNU.f
  )
target_link_libraries(${PROJECT_NAME} LBT_LIB)

# install compiled binary
install (TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_NAME})
