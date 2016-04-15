cmake_minimum_required(VERSION 2.8)

project(Vpr)

#set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

#include(cotire)
message(STATUS ${CMAKE_CURRENT_BINARY_DIR})
message(STATUS ${CMAKE_BINARY_DIR})

add_subdirectory(../libarchfpga ../libarchfpga)
add_subdirectory(../pcre ../pcre)
add_subdirectory(../printhandler ../printhandler)

set(srcs
	SRC/base/CheckArch.c
	SRC/base/check_netlist.c
	SRC/base/CheckOptions.c
	SRC/base/CheckSetup.c
	SRC/base/draw.c
	SRC/base/globals.c
	SRC/base/graphics.c
	SRC/base/OptionTokens.c
	SRC/base/place_and_route.c
	SRC/base/read_blif.c
	SRC/base/read_netlist.c
	SRC/base/ReadOptions.c
	SRC/base/read_place.c
	SRC/base/read_settings.c
	SRC/base/SetupGrid.c
	SRC/base/SetupVPR.c
	SRC/base/ShowSetup.c
	SRC/base/stats.c
	SRC/base/verilog_writer.c
	SRC/base/vpr_api.c
	SRC/base/barrier.c
	SRC/base/bounding_box.c
	SRC/pack/cluster.c
	SRC/pack/cluster_feasibility_filter.c
	SRC/pack/cluster_legality.c
	SRC/pack/cluster_placement.c
	SRC/pack/output_blif.c
	SRC/pack/output_clustering.c
	SRC/pack/pack.c
	SRC/pack/pb_type_graph_annotations.c
	SRC/pack/pb_type_graph.c
	SRC/pack/prepack.c
	SRC/pack/print_netlist.c
	SRC/place/place.c
	SRC/place/place_macro.c
	SRC/place/place_stats.c
	SRC/place/timing_place.c
	SRC/place/timing_place_lookup.c
	SRC/power/power.c
	SRC/power/power_callibrate.c
	SRC/power/power_cmos_tech.c
	SRC/power/power_components.c
	SRC/power/power_lowlevel.c
	SRC/power/power_sizing.c
	SRC/power/PowerSpicedComponent.c
	SRC/power/power_util.c
	SRC/route/check_route.c
	SRC/route/check_rr_graph.c
	SRC/route/route_breadth_first.c
	SRC/route/route_common.c
	SRC/route/route_timing.c
	SRC/route/route_tree_timing.c
	SRC/route/rr_graph2.c
	SRC/route/rr_graph_area.c
	SRC/route/rr_graph.c
	SRC/route/rr_graph_indexed_data.c
	SRC/route/rr_graph_sbox.c
	SRC/route/rr_graph_timing_params.c
	SRC/route/rr_graph_util.c
	SRC/route/segment_stats.c
	#SRC/route/parallel_route_timing.c
	#SRC/route/new_cost_parallel_route_timing.c
	#SRC/route/fine_grained_parallel_route_timing.c
	#SRC/route/fine_grained_parallel_route_timing_2.c
	#SRC/route/test.c
	#SRC/route/partition.cxx
	SRC/parallel_route/route_tree_timing_local.c
	SRC/parallel_route/route_common_local.c
	SRC/parallel_route/fine.cxx
	#SRC/parallel_route/phased.cxx
	SRC/parallel_route/congestion.cxx
	SRC/parallel_route/router.cxx
	SRC/parallel_route/init.cxx
	SRC/parallel_route/log.cxx
	SRC/parallel_route/spatial.cxx
	SRC/parallel_route/astar.cxx
	#SRC/parallel_route/partitioning_route.c
	SRC/parallel_route/route_tree.c
	#SRC/parallel_route/trace.c
	SRC/parallel_route/utility.c
	SRC/parallel_route/net_cluster.cxx
	SRC/parallel_route/scheduler.c
	#SRC/parallel_route/advanced_parallel_route_timing.c
	SRC/timing/net_delay.c
	SRC/timing/path_delay2.c
	SRC/timing/path_delay.c
	#SRC/timing/path_delay_new.c
	SRC/timing/read_sdc.c
	SRC/timing/slre.c
	SRC/util/hash.c
	SRC/util/heapsort.c
	SRC/util/token.c
	SRC/util/vpr_utils.c
	)

set_source_files_properties(${srcs} PROPERTIES LANGUAGE CXX)

add_library(Vpr STATIC
	${srcs}
)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	set(external_includes 
		/usr/include/libxml2
		/usr/local/include
		/opt/intel/vtune_amplifier_xe/include/
		/opt/intel/tbb/include/
		/home/chinhau5/Downloads/boost_1_60_0/
		/home/chinhau5/Downloads/tbb44_20151115oss/include
		/home/chinhau5/Downloads/mlpack-1.0.12/build/include
		/home/chinhau5/Downloads/armadillo-6.400.3/include
		/home/chinhau5/Downloads/metis-5.1.0/include/
		)
else()
	set(external_includes 
		/usr/local/include
		/opt/local/include
		/usr/include/libxml2
		#/usr/local/Cellar/mpich/3.2/include
		/usr/local/Cellar/open-mpi/1.10.2/include
		/opt/X11/include
		/Users/chinhau5/Downloads/tbb44_20151115oss/include
		/Users/chinhau5/Downloads/mlpack-1.0.12/build/include
		/Users/chinhau5/Downloads/armadillo-6.400.3/include
		/Users/chinhau5/Downloads/metis-5.1.0/include/
		)
	#set(pch_flags
		#-fasm-blocks
		#-fretain-comments-from-system-headers
		#)
endif()

target_include_directories(Vpr PUBLIC 
	${external_includes}
	SRC/util
	SRC/timing 
	SRC/pack 
	SRC/place 
	SRC/base 
	SRC/route 
	SRC/power 
	SRC/parallel_route
)

target_compile_options(Vpr PRIVATE
	-Wno-write-strings -Wno-unused-variable -Wno-unused-private-field
	#-DNO_GRAPHICS
	) 

target_link_libraries(Vpr
	ArchFpga
	PrintHandler
	Pcre
	)

if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
else()
	#get_property(dirs TARGET Vpr PROPERTY INCLUDE_DIRECTORIES)
	#foreach(dir ${dirs})
		#set(include_dirs ${include_dirs} -I ${dir})
	#endforeach()

	#get_property(opts TARGET Vpr PROPERTY COMPILE_OPTIONS)
	#get_property(flags TARGET Vpr PROPERTY COMPILE_FLAGS)

	#message(STATUS ${dirs})
	#message(STATUS ${flags})

	#add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pch.h.gch
	#COMMAND ${CMAKE_CXX_COMPILER} -x c++-header ${pch_flags} ${opts} ${include_dirs} ${CMAKE_CURRENT_SOURCE_DIR}/SRC/route/pch.h -o ${CMAKE_CURRENT_BINARY_DIR}/pch.h.gch
	#DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/SRC/route/pch.h)

	#add_custom_target(Pch
	#DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pch.h.gch)

	#add_dependencies(Vpr Pch)
	#target_compile_options(Vpr PRIVATE
	#-include-pch ${CMAKE_CURRENT_BINARY_DIR}/pch.h.gch
	#)
endif()

#cotire(Vpr)
