cmake_minimum_required(VERSION 3.10)
project(Muography LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_library(MuographyLib INTERFACE)
target_include_directories(MuographyLib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(MuographyLib INTERFACE ${ROOT_LIBRARIES} ROOT::Minuit)

add_executable(SimulationWaterPipe work/SimulationWaterPipe.cxx)
target_link_libraries(SimulationWaterPipe MuographyLib)

add_executable(ReconstructionWaterPipe work/ReconstructionWaterPipe.cxx)
target_link_libraries(ReconstructionWaterPipe MuographyLib)

add_executable(SimulationFiberDetector work/SimulationFiberDetector.cxx)
target_link_libraries(SimulationFiberDetector MuographyLib)

add_executable(ReconstructionFiberDetector work/ReconstructionFiberDetector.cxx)
target_link_libraries(ReconstructionFiberDetector MuographyLib)

add_executable(SimulationKomen work/SimulationKomen.cxx)
target_link_libraries(SimulationKomen MuographyLib)

add_executable(ReconstructionKomen work/ReconstructionKomen.cxx)
target_link_libraries(ReconstructionKomen MuographyLib)

add_custom_target(tidy
    COMMAND rm -f ${CMAKE_CURRENT_SOURCE_DIR}/*~
    COMMAND rm -f ${CMAKE_CURRENT_SOURCE_DIR}/*/*~
    COMMAND rm -f ${CMAKE_CURRENT_SOURCE_DIR}/*/*/*~
    COMMAND rm -rf ${CMAKE_BINARY_DIR}/bin/*
    COMMENT "Cleaning up backup files and binaries..."
)
