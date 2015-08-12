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
	-Wno-write-strings -Wno-unused-variable -Wno-unused-private-field)
