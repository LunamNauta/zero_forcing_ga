cmake_minimum_required(VERSION 3.15)
project(ZeroForcing VERSION 1.0 LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT GUROBI_DIR)
  set(GUROBI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../gurobi/gurobi1103/linux64")
endif()
if(NOT GUROBI_INC)
    set(GUROBI_INC "${GUROBI_DIR}/include")
endif()
if(NOT GUROBI_LIB_DIR)
    set(GUROBI_LIB_DIR "${GUROBI_DIR}/lib")
endif()
message(STATUS "Gurobi Directory: ${GUROBI_DIR}")
message(STATUS "Gurobi Include:   ${GUROBI_INC}")
message(STATUS "Gurobi Lib Dir:   ${GUROBI_LIB_DIR}")

add_compile_options(-O3 -march=native -mtune=native -lto)

include_directories(include)
include_directories(${GUROBI_INC})

file(GLOB_RECURSE LIB_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(graph_lib STATIC ${LIB_SOURCES})

file(GLOB EXAMPLE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/examples/*.cpp")

foreach(FILE_PATH ${EXAMPLE_FILES})
  get_filename_component(BIN_NAME ${FILE_PATH} NAME_WE)
  add_executable(${BIN_NAME} ${FILE_PATH})

  target_link_libraries(${BIN_NAME} PRIVATE 
    graph_lib                   
    "${GUROBI_LIB_DIR}/libgurobi_c++.a"
    "${GUROBI_LIB_DIR}/libgurobi110.so"
    m 
    dl 
    pthread
    nauty
  )

  set_target_properties(${BIN_NAME} PROPERTIES 
    INSTALL_RPATH "${GUROBI_LIB_DIR}"
    BUILD_WITH_INSTALL_RPATH TRUE
  )
endforeach()