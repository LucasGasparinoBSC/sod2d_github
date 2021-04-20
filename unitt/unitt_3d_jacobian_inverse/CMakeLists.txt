cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

get_filename_component(projectId ${CMAKE_CURRENT_LIST_DIR} NAME)
#string(REPLACE " " "_" ProjectId ${ProjectId})

project(${projectId} LANGUAGES Fortran)
#file(GLOB_RECURSE SRC_FILES ${CMAKE_HOME_DIRECTORY}/src/*.f90)
file(GLOB_RECURSE SRC_FILES ${CMAKE_HOME_DIRECTORY}/src/jacobian_oper.f90 ${CMAKE_HOME_DIRECTORY}/unitt/unitt_3d_jacobian_inverse/*.f90)

add_executable(${PROJECT_NAME} ${SRC_FILES})   # Build executable from listed sources

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
