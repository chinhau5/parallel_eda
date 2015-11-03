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
	-ansi -pedantic
	-Wno-write-strings -Wno-unused-variable -Wno-unused-private-field
	-DEZXML_NOMMAP -D_POSIX_C_SOURCE)
