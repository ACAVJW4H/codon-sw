cmake_minimum_required (VERSION 2.8.1)

project(CodonAligner)

set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -pedantic")
include_directories("deps/seqan/core/include")
include_directories("deps/seqan/extras/include")
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#list(APPEND CMAKE_MODULE_PATH 
  #${CMAKE_CURRENT_SOURCE_DIR}/deps/seqan/util/cmake)
#find_package(SeqAn)


set(CPP_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/seqan_codons.cpp
  )
add_library(codonaligner-static STATIC ${CPP_FILES})

add_executable(codon-align ${CMAKE_CURRENT_SOURCE_DIR}/src/codon_align.cpp)

find_package(ZLIB REQUIRED)
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( codon-align ${ZLIB_LIBRARIES} )
    set(SEQAN_HAS_ZLIB true)
    add_definitions( -DSEQAN_HAS_ZLIB)
endif( ZLIB_FOUND )

find_package(BZip2 REQUIRED)
if( BZIP2_FOUND )
    include_directories( ${BZIP2_INCLUDE_DIRS} )
    target_link_libraries( codon-align ${BZIP2_LIBRARIES} )
    set(SEQAN_HAS_BZIP2 true)
    add_definitions( -DSEQAN_HAS_BZIP2 )

endif( BZIP2_FOUND )

target_link_libraries(codon-align codonaligner-static)
install(TARGETS codon-align
        RUNTIME DESTINATION bin)
