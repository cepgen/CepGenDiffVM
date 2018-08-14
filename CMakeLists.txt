cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(CepGenDiffVM)

set(CEPGEN_PATH ${CEPGEN_PATH} $ENV{CEPGEN_PATH})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${CEPGEN_PATH}/cmake)
find_package(CepGen REQUIRED)

file(GLOB sources src/*.cpp)

add_library(CepGenDiffVM SHARED ${sources})
target_link_libraries(CepGenDiffVM PRIVATE CepGen::CepGen)
target_include_directories(CepGenDiffVM PRIVATE ${PROJECT_SOURCE_DIR})

# copy the input cards
file(GLOB_RECURSE input_cards RELATIVE ${PROJECT_SOURCE_DIR} cards/*)
foreach(_files ${input_cards})
  configure_file(${_files} ${_files} COPYONLY)
endforeach()

