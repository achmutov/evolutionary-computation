cmake_minimum_required(VERSION 3.15)
project(evolutionary_computation VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(CMAKE_BUILD_TYPE MATCHES Release)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(
            -O3
            -march=native
            -ffast-math
            -fno-exceptions
            -fno-rtti
            -DNDEBUG
        )
        add_link_options(-flto)
    endif()

    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported)
    if(ipo_supported)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
else()
    add_compile_definitions(EXCEPTIONS=1)
endif()

include_directories(include)
add_subdirectory(src)
add_subdirectory(samples)
