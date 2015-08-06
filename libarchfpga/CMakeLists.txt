cmake_minimum_required(VERSION 2.8)

project(ArchFpga)

set(srcs
	read_xml_arch_file.c 
	read_xml_util.c 
	ezxml.c 
	ReadLine.c 
	util.c
	)

set_source_files_properties(${srcs} PROPERTIES LANGUAGE CXX)

add_library(ArchFpga STATIC
	${srcs}
)

target_include_directories(ArchFpga PUBLIC
	include
	../printhandler/SRC/TIO_InputOutputHandlers
)

target_compile_options(ArchFpga PRIVATE
	-std=c++11
	-g -O0
	-fno-omit-frame-pointer #-fPIC
	-Wall -Wno-write-strings -Wno-unused-variable -Wno-unused-private-field)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND SANITIZE EQUAL 1)
target_compile_options(ArchFpga PRIVATE
	-fsanitize=memory
	-fsanitize-blacklist=/home/chinhau5/vpr_parallel/blacklist.txt)
endif()
