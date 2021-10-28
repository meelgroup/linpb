cmake_minimum_required(VERSION 3.5)
project(LinPB)
message(STATUS "########## LinPB ###########")
message(STATUS "build type is ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 17)

# Feel free to decrease the required boost version if you can verify that it works with
# a smaller version. It is known not to compile with version 1.53.0 and lower.
find_package(Boost 1.67 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")

if (NOT CMAKE_BUILD_TYPE)
    message(WARNING "Please use -DCMAKE_BUILD_TYPE=[build_type] where [build_type] is one of: Debug, Release, RelWithDebInfo and MinSizeRel")
endif()

set(build_static "OFF"
    CACHE STRING "Build and use static libraries.")

set(profile_flags "OFF"
    CACHE STRING "Build using gprof flag (-pg).")

set(soplex "OFF"
    CACHE STRING "Use SoPlex to check for rational infeasibility.")

set(soplex_pkg "${PROJECT_SOURCE_DIR}/soplex-5.0.2.tgz"
    CACHE STRING "SoPlex package, download at https://soplex.zib.de/download.php?fname=soplex-5.0.2.tgz")

set(gmp "OFF"
        CACHE STRING "Use GMP instead of Boost's builtin multiple precision integer.")

set(testruns_timeout 1
    CACHE STRING "Timeout of individual test runs.")

IF(soplex)
    # provides libsoplex and ${SOPLEX_INCLUDE_DIRS}
    include(${PROJECT_SOURCE_DIR}/cmake/soplex_build_and_load.cmake)
endif()

if(gmp)
    find_package(GMP REQUIRED)
    include_directories(${GMP_INCLUDE_DIR})
endif()

# Look for M4RI used by Gauss Jordan Elimination
option(REQUIRE_M4RI "Must use m4ri" OFF)
option(NOM4RI "Don't use m4ri" OFF)

if (NOT NOM4RI)
    find_package(M4RI)
    IF (M4RI_FOUND)
        MESSAGE(STATUS "OK, Found M4RI lib at ${M4RI_LIBRARIES} and includes at ${M4RI_INCLUDE_DIRS}")
        add_definitions( -DUSE_M4RI )
    ELSE (M4RI_FOUND)
        MESSAGE(WARNING "Did not find M4RI, XOR detection&manipulation disabled")
        if (REQUIRE_M4RI)
            MESSAGE(FATAL_ERROR "REQUIRE_M4RI was set but M4RI was not found!")
        endif()
    ENDIF (M4RI_FOUND)
endif()

if(${build_static})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -static")
endif()

if(${profile_flags})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

# Get the current working branch
execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# Get the latest abbreviated commit hash of the working branch
execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")

# add -DUSE_GAUSS
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++1z -DUSE_GAUSS")

set(source_files
    src/Constr.cpp
    src/ConstrExp.cpp
    src/ConstrSimple.cpp
    src/LpSolver.cpp
    src/Solver.cpp
    src/SolverStructs.cpp
    src/parsing.cpp
    src/quit.cpp
    src/roundingsat.cpp
    src/run.cpp
    # add gauss files
    src/gaussian/solverwrapper.cpp
    src/gaussian/gaussian.cpp
    src/gaussian/packedrow.cpp
)

set(header_files
    src/Constr.hpp
    src/ConstrExp.hpp
    src/ConstrSimple.hpp
    src/SolverStructs.hpp
    src/Logger.hpp
    src/Stats.hpp
    src/LpSolver.hpp
    src/Options.hpp
    src/IntSet.hpp
    src/Solver.hpp
    src/typedefs.hpp
    src/aux.hpp
    src/globals.hpp
    src/parsing.hpp
    src/quit.hpp
    src/run.hpp
    # add gauss files
    src/gaussian/Vec.h
    src/gaussian/constants.h
    src/gaussian/gausswatched.h
    src/gaussian/propby.h
    src/gaussian/xor.h
    src/gaussian/XAlloc.h
    src/gaussian/gqueuedata.h
    src/gaussian/packedrow.h
    src/gaussian/solvertypes.h
    src/gaussian/solverwrapper.h
    src/gaussian/bitarray.h
    src/gaussian/gaussian.h
    src/gaussian/packedmatrix.h
    src/gaussian/popcnt.h
    src/gaussian/solvertypesmini.h.in
    src/gaussian/vardata.h
    src/gaussian/toplevelgaussabst.h
)

# m4ri
if (M4RI_FOUND)
	set(source_files ${source_files} src/gaussian/toplevelgauss.cpp)
	set(header_files ${header_files} src/gaussian/toplevelgauss.h)
endif (M4RI_FOUND)

set(all_files ${source_files} ${header_files})

# add_executable(linpb ${source_files})
##################################################
# add linpb as a library
add_library(linpb STATIC ${source_files})
set(PROJECT_VERSION_MAJOR 2)
set(PROJECT_VERSION_MINOR 0)
set_target_properties(linpb PROPERTIES
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)
install(TARGETS linpb
    EXPORT ${CRYPTOMINISAT5_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/approxmc"
)#################################################

target_include_directories(linpb PRIVATE ${PROJECT_SOURCE_DIR}/src/ ${SOPLEX_INCLUDE_DIRS})

IF(soplex)
    target_link_libraries(linpb libsoplex)
    target_compile_definitions(linpb PUBLIC WITHSOPLEX)
endif()

IF(gmp)
    target_link_libraries(linpb gmp gmpxx)
    target_compile_definitions(linpb PUBLIC WITHGMP)
endif()

if (M4RI_FOUND)
    include_directories(${M4RI_INCLUDE_DIRS})
    target_link_libraries(linpb ${M4RI_LIBRARIES})
endif (M4RI_FOUND)

FIND_PROGRAM(CLANG_FORMAT "clang-format")
IF(CLANG_FORMAT)
    add_custom_target(
        format
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMAND ${CLANG_FORMAT}
        -style=file
        -i
        ${all_files}
    )
endif()

add_custom_target(
    testruns
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
    COMMAND ./run_tests.sh ${testruns_timeout} "testruns" ${CMAKE_CURRENT_BINARY_DIR}/linpb
)
add_dependencies(testruns linpb)
message(STATUS "########### end ############")
