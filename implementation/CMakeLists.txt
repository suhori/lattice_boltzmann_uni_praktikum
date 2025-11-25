cmake_minimum_required(VERSION 4.0)
project(lattice_boltzmann_uni_praktikum)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(.)

add_executable(lattice_boltzmann_uni_praktikum
        LBM.cpp
        LBM.h
        main.cpp
        seconds.cpp
        seconds.h)

# If building with Clang, prefer libc++ over libstdc++ so that C++23 features like std::mdspan are available.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Use libc++ standard library implementation
    target_compile_options(lattice_boltzmann_uni_praktikum PRIVATE -stdlib=libc++)
    target_link_options(lattice_boltzmann_uni_praktikum PRIVATE -stdlib=libc++)

    # On Linux, explicitly link against libc++ and libc++abi to avoid picking up libstdc++ implicitly.
    if (UNIX AND NOT APPLE)
        target_link_libraries(lattice_boltzmann_uni_praktikum PRIVATE c++ c++abi)
    endif ()
endif ()
