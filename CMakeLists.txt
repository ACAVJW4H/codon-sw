cmake_minimum_required (VERSION 2.8.1)

project(CodonAligner)

# Versioning
set(CMAKE_MODULE_PATH
    ${CMAKE_CURRENT_SOURCE_DIR}/deps/cmake
    ${CMAKE_MODULE_PATH})
include(GetGitRevisionDescription)

git_describe(GIT_DESCRIPTION --always --long)
configure_file (
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )
include_directories("${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -Wextra -pedantic")
set(SEQAN_C++11_STANDARD "ON")

# Flags for Debug / Release
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  add_definitions(-DSEQAN_ENABLE_DEBUG=0)
endif()
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_definitions(-DSEQAN_ENABLE_DEBUG=1)
endif()

include_directories("deps/seqan/core/include")
include_directories("deps/seqan/extras/include")
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CPP_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/seqan_codons.cpp
  )
add_library(codonaligner-static STATIC ${CPP_FILES})

add_executable(codon-sw ${CMAKE_CURRENT_SOURCE_DIR}/src/codon_align.cpp)

# ZLIB required for BAM
find_package(ZLIB REQUIRED)
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( codon-sw ${ZLIB_LIBRARIES} )
    set(SEQAN_HAS_ZLIB true)
    add_definitions(-DSEQAN_HAS_ZLIB)
endif( ZLIB_FOUND )

# Optional: for BZIP-compressed sequences
find_package(BZip2)
if( BZIP2_FOUND )
    include_directories( ${BZIP2_INCLUDE_DIRS} )
    target_link_libraries( codon-sw ${BZIP2_LIBRARIES} )
    set(SEQAN_HAS_BZIP2 true)
    add_definitions( -DSEQAN_HAS_BZIP2 )

endif( BZIP2_FOUND )

target_link_libraries(codon-sw codonaligner-static)
install(TARGETS codon-sw
        RUNTIME DESTINATION bin)
