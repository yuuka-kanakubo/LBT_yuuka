set_source_files_properties(JetOutput-0.7.f PROPERTIES COMPILE_FLAGS "-w")
set_source_files_properties(hydro-CCNU.f PROPERTIES COMPILE_FLAGS "-w")
add_library(
		LBT_LIB STATIC
		main_functions.cpp
		LBTConfig.cpp
		LBTcl.cpp
		JetOutput-0.7.f
		hydro-CCNU.f
	   )

set(
    CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_HOME_DIRECTORY}"
   )
add_executable("${PROJECT_NAME}"
  main.cpp
  )
target_link_libraries("${PROJECT_NAME}" LBT_LIB)


install (TARGETS "${PROJECT_NAME}" DESTINATION "${PROJECT_NAME}")
