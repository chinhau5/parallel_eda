cmake_minimum_required(VERSION 2.8)

project(Router)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	link_directories(
		/usr/local/lib
		/usr/X11/lib
		/opt/intel/tbb/lib/intel64/gcc4.4
		/home/chinhau5/Downloads/metis-5.1.0/build/Linux-x86_64/libmetis
		)
else()
	link_directories(
		/usr/local/lib
		/opt/local/lib
		/usr/local/Cellar/mpich/3.2/lib
		/Users/chinhau5/Downloads/tbb44_20151115oss/build/macos_intel64_clang_cc5.2.0_os10.11.1_debug
		/Users/chinhau5/Downloads/tbb44_20151115oss/build/macos_intel64_clang_cc5.2.0_os10.11.1_release
		/Users/chinhau5/Downloads/mlpack-1.0.12/build/lib
		/Users/chinhau5/Downloads/metis-5.1.0/build/Darwin-x86_64/libmetis/Debug
		)
endif()

#set(CMAKE_SKIP_RPATH true)
#-fprofile-arcs -ftest-coverage -fsanitize=address 
#set(CMAKE_CXX_FLAGS "-std=c++11 -g -O0 -pthread -fno-omit-frame-pointer -Wall -Wno-write-strings -Wno-unused-variable -Wno-unused-private-field -DNO_GRAPHICS") 

#generic compile options go here
add_compile_options(
	#-D TBB_USE_DEBUG=1
	#-ferror-limit=0
	-ftemplate-backtrace-limit=0
	-std=c++14
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

add_subdirectory(vpr)

set_source_files_properties(vpr/SRC/main.c PROPERTIES LANGUAGE CXX)

add_executable(Router vpr/SRC/main.c)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	set(link_libs
		Vpr
		X11
		pthread
		dl 
		m
		zlog
		boost_timer
		boost_system
		tbb
		metis
		#-Wl,--eh-frame-hdr
		)
else()
	set(link_libs
		Vpr
		mpicxx 
		mpi 
		pmpi
		X11
		pthread
		dl 
		m
		zlog
		boost_timer-mt-d
		boost_system-mt-d
		igraph
		tbb
		mlpack
		metis
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

set_target_properties(Router PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
