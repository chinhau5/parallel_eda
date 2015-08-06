cmake_minimum_required(VERSION 2.8)

project(Pcre)

set(srcs
	SRC/get.c
	SRC/pcre.c
	)

set_source_files_properties(${srcs} PROPERTIES LANGUAGE CXX)

add_library(Pcre STATIC
	${srcs}
	)

target_compile_options(Pcre PRIVATE
	-std=c++11
	-g -O0
	-ansi -pedantic
	-fno-omit-frame-pointer #-fPIC
	-Wall -Wno-write-strings -Wno-unused-variable -Wno-unused-private-field
	-D EZXML_NOMMAP -D_POSIX_C_SOURCE)

#target_include_directories(ArchFpga PUBLIC
#	${CMAKE_CURRENT_SOURCE_DIR}/include
#	../printhandler/SRC/TIO_InputOutputHandlers 
#	)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND SANITIZE EQUAL 1)
target_compile_options(Pcre PRIVATE
	-fsanitize=memory
	-fsanitize-blacklist=/home/chinhau5/vpr_parallel/blacklist.txt)
endif()

