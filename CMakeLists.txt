cmake_minimum_required(VERSION 2.8)

project(Router)

link_directories(
	/usr/local/lib
	/opt/local/lib
	~/Downloads/tbb44_20151115oss/build/macos_intel64_clang_cc5.2.0_os10.11.1_debug
	~/Downloads/tbb44_20151115oss/build/macos_intel64_clang_cc5.2.0_os10.11.1_release
	/Users/chinhau5/Downloads/mlpack-1.0.12/build/lib
	)

#set(CMAKE_SKIP_RPATH true)
#-fprofile-arcs -ftest-coverage -fsanitize=address 
#set(CMAKE_CXX_FLAGS "-std=c++11 -g -O0 -pthread -fno-omit-frame-pointer -Wall -Wno-write-strings -Wno-unused-variable -Wno-unused-private-field -DNO_GRAPHICS") 

#generic compile options go here
add_compile_options(
	#-D TBB_USE_DEBUG=1
	-ferror-limit=0
	-std=c++11
	-g -O0
	#-Bdynamic
	-Wall
	#	-H
	)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND SANITIZE EQUAL 1)
add_compile_options(
	-fsanitize=thread
	#-fsanitize-blacklist=/home/chinhau5/vpr_parallel/blacklist.txt
	)
endif()

add_subdirectory(libarchfpga)
add_subdirectory(pcre)
add_subdirectory(printhandler)
add_subdirectory(vpr)

set_source_files_properties(vpr/SRC/main.c PROPERTIES LANGUAGE CXX)

add_executable(Router vpr/SRC/main.c)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	set(link_libs
		Vpr
		pthread
		dl 
		m
		zlog
		boost_timer
		boost_system
		tbb_debug
		)
else()
	set(link_libs
		Vpr
		pthread
		dl 
		m
		zlog
		boost_timer-mt-d
		boost_system-mt-d
		igraph
		tbb
		mlpack
		)
endif()

target_link_libraries(Router 
	${link_libs}
	)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND SANITIZE EQUAL 1)
target_link_libraries(Router
	-fsanitize=thread
	#-fsanitize-blacklist=/home/chinhau5/vpr_parallel/blacklist.txt
	)
endif()

set_target_properties(Router PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
