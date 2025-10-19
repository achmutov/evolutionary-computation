add_library(evolutionary_computation)

target_link_directories(evolutionary_computation PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>)

target_sources(evolutionary_computation PRIVATE cli.cpp)
add_subdirectory(loader)
add_subdirectory(solver)
