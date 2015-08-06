cmake_minimum_required(VERSION 2.8)

project(Router)

#link_directories(
#	${CMAKE_CURRENT_SOURCE_DIR}/libarchfpga
#	)

set(CMAKE_SKIP_RPATH true)
#-fprofile-arcs -ftest-coverage -fsanitize=address 
#set(CMAKE_CXX_FLAGS "-std=c++11 -g -O0 -pthread -fno-omit-frame-pointer -Wall -Wno-write-strings -Wno-unused-variable -Wno-unused-private-field -DNO_GRAPHICS") 

add_subdirectory(libarchfpga)
add_subdirectory(pcre)
add_subdirectory(printhandler)
add_subdirectory(vpr)

set_source_files_properties(vpr/SRC/main.c PROPERTIES LANGUAGE CXX)

add_executable(Router vpr/SRC/main.c)

target_link_libraries(Router 
	Vpr
	pthread
	dl 
	m
	)

target_compile_options(Router PRIVATE
	-std=c++11
	-g -O0
	-fno-omit-frame-pointer
	-Wall)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND SANITIZE EQUAL 1)
target_compile_options(Router PRIVATE
	-fsanitize=memory
	-fsanitize-blacklist=/home/chinhau5/vpr_parallel/blacklist.txt)

target_link_libraries(Router
	-fsanitize=memory
	-fsanitize-blacklist=/home/chinhau5/vpr_parallel/blacklist.txt)
endif()

set_target_properties(Router PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
