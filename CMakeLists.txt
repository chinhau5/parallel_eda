cmake_minimum_required(VERSION 2.8)

project(Router)

option(WITH_NUMA "Compile and link with NUMA library" OFF)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	link_directories(
		/usr/local/lib
		/usr/X11/lib
		#/opt/intel/tbb/lib/intel64/gcc4.4
		/home/chinhau5/Downloads/metis-5.1.0/build/Linux-x86_64/libmetis
		/home/chinhau5/Downloads/tbb44_20160526oss/lib/intel64/gcc4.4
		/home/chinhau5/zlog/lib
		/home/chinhau5/libs/lib
		#/home/chinhau/Downloads/metis-5.1.0/build/Linux-x86_64/libmetis
		#/home/chinhau/Downloads/tbb44_20160526oss/lib/intel64/gcc4.4/
		#/home/chinhau/zlog/lib
		/home/akkumar/Downloads/metis-5.1.0/build/Linux-x86_64/libmetis
		/home/akkumar/Downloads/tbb44_20160526oss/lib/intel64/gcc4.4
		/home/akkumar/Downloads/boost_1_60_0/stage/lib
		/home/akkumar/zlog/lib
		/home/akkumar/libs/lib
		#/home/svu/a0045424/Downloads/metis-5.1.0/build/Linux-x86_64/libmetis
		#/home/svu/a0045424/Downloads/tbb44_20160526oss/lib/intel64/gcc4.4
		#/home/svu/a0045424/Downloads/boost_1_60_0/stage/lib
		#/home/svu/a0045424/zlog/lib
		)
else()
	link_directories(
		/usr/local/lib
		/opt/X11/lib
		/Users/chinhau5/Downloads/tbb44_20151115oss/build/macos_intel64_clang_cc5.2.0_os10.11.1_debug
		/Users/chinhau5/Downloads/tbb44_20151115oss/build/macos_intel64_clang_cc5.2.0_os10.11.1_release
		/Users/chinhau5/Downloads/metis-5.1.0/build/Darwin-x86_64/libmetis/
		/Users/chinhau5/libs/lib
		)
endif()

#set(CMAKE_SKIP_RPATH true)
#-fprofile-arcs -ftest-coverage -fsanitize=address 
#set(CMAKE_CXX_FLAGS "-std=c++11 -g -O0 -pthread -fno-omit-frame-pointer -Wall -Wno-write-strings -Wno-unused-variable -Wno-unused-private-field -DNO_GRAPHICS") 

#generic compile options go here
add_compile_options(
	#-DTBB_USE_DEBUG=1
	#-ferror-limit=0
	-ftemplate-backtrace-limit=0
	-std=c++11
	-g -O3
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
		#tbbmalloc_proxy
		X11
		pthread
		dl 
		m
		rt
		zlog
		boost_timer
		boost_system
		boost_thread
		tbb
		metis
		#mlpack
		)
else()
	set(link_libs
		Vpr
		mpi 
		X11
		pthread
		dl 
		m
		zlog
		boost_timer-mt
		boost_system-mt
		boost_thread-mt
		tbb
		metis
		mlpack
		)
endif()

target_link_libraries(Router 
	${link_libs}
	)

if (WITH_NUMA)
	target_link_libraries(Router 
		numa
		)
endif()

if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND SANITIZE EQUAL 1)
target_link_libraries(Router
	-fsanitize=thread
	#-fsanitize-blacklist=/home/chinhau5/vpr_parallel/blacklist.txt
	)
endif()

set_target_properties(Router PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
