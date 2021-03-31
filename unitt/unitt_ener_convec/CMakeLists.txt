cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

get_filename_component(sod2d ${CMAKE_CURRENT_LIST_DIR} NAME)

project(sod2d LANGUAGES Fortran)
#file(GLOB_RECURSE SRC_FILES ${CMAKE_HOME_DIRECTORY}/src/*.f90)
file(GLOB_RECURSE SRC_FILES ${CMAKE_HOME_DIRECTORY}/src/elem_convec.f90 ${CMAKE_HOME_DIRECTORY}/unitt/unitt_ener_convec/*.f90)

add_executable(unitt.x ${SRC_FILES})   # Build executable from listed sources

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
