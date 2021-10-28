if(NOT EXISTS ${soplex_pkg})
    message(FATAL_ERROR "Please provide location of SoPlex package via -Dsoplex_pkg=..., did not find it at: ${soplex_pkg}.")
endif()

set(soplex_dir "${PROJECT_BINARY_DIR}/soplex")
set(soplex_targets "${soplex_dir}/install/lib/cmake/soplex/soplex-targets.cmake")
set(soplex_config "${soplex_dir}/install/lib/cmake/soplex/soplex-config.cmake")

# Unpack and compile soplex at configure time
if(NOT EXISTS ${soplex_targets} OR NOT EXISTS ${soplex_config})
    if(${build_static})
        set(soplex_cmake_args
            ${soplex_cmake_args}
            # set CMAKE_FIND_LIBRARY_SUFFIXES ".a"
            -DCMAKE_PROJECT_SOPLEX_INCLUDE=${PROJECT_SOURCE_DIR}/cmake/soplex_Settings_Overwrite.cmake
            -DCMAKE_EXE_LINKER_FLAGS=-static-libgcc -static-libstdc++ -static
            -DBUILD_SHARED_LIBS=OFF
            -DBoost_USE_STATIC_LIBS=ON
        )
    endif()

    set(soplex_cmake_args
        ${soplex_cmake_args}
        -DCMAKE_BUILD_TYPE=Release
        -DCMAKE_INSTALL_PREFIX=${soplex_dir}/install
        -DBOOST=ON
        -DGMP=OFF
        -DMPFR=OFF
        -DBUILD_TESTING=OFF
        -DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=OFF
        -DZLIB=OFF
    )

    # TODO: linking with SoPlex seems to use the build type of RoundingSat, no matter the build type we configure above?

    configure_file(${PROJECT_SOURCE_DIR}/cmake/soplex_CMakeLists.txt.in ${soplex_dir}/download/CMakeLists.txt)

    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${soplex_dir}/download)
    if(result)
      message(FATAL_ERROR "CMake step for SoPlex failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build . -j
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${soplex_dir}/download)
    if(result)
      message(FATAL_ERROR "Build step for SoPlex failed: ${result}")
    endif()
endif()

# provides libsoplex
include(${soplex_dir}/install/lib/cmake/soplex/soplex-targets.cmake)
# provides ${SOPLEX_INCLUDE_DIRS}
include(${soplex_dir}/install/lib/cmake/soplex/soplex-config.cmake)
