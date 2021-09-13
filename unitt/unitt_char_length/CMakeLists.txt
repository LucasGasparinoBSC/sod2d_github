cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

get_filename_component(projectId ${CMAKE_CURRENT_LIST_DIR} NAME)
#string(REPLACE " " "_" ProjectId ${ProjectId})

project(${projectId} LANGUAGES Fortran)
file(GLOB SRC_FILES ${CMAKE_HOME_DIRECTORY}/src/mod_geom.f90 ${CMAKE_HOME_DIRECTORY}/src/elem_qua.f90 ${CMAKE_HOME_DIRECTORY}/src/elem_hex.f90 ${CMAKE_HOME_DIRECTORY}/unitt/unitt_char_length/*.f90)

add_executable(${PROJECT_NAME} ${SRC_FILES})   # Build executable from listed sources

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
