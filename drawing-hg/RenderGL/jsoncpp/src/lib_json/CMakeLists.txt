set(CMAKE_VERBOSE_MAKEFILE on)

# the variable "libjsoncpp_SRCS" contains all .cpp files of this project
set(libjsoncpp_SRCS
    json_reader.cpp
    json_value.cpp
    json_writer.cpp
    json_batchallocator.h
)

add_library(jsoncpp STATIC ${libjsoncpp_SRCS})
include_directories(jsoncpp ${} ${CMAKE_SOURCE_DIR}/include )
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/../lib)
