cmake_minimum_required(VERSION 2.8)

project(Vpr)

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
	SRC/route/route_common_local.c
	SRC/route/route_timing.c
	SRC/route/parallel_route_timing.c
	SRC/route/new_cost_parallel_route_timing.c
	SRC/route/fine_grained_parallel_route_timing.c
	SRC/route/fine_grained_parallel_route_timing_2.c
	SRC/route/hybrid_route.c
	SRC/route/advanced_parallel_route_timing.c
	SRC/route/route_tree_timing.c
	SRC/route/route_tree_timing_local.c
	SRC/route/rr_graph2.c
	SRC/route/rr_graph_area.c
	SRC/route/rr_graph.c
	SRC/route/rr_graph_indexed_data.c
	SRC/route/rr_graph_sbox.c
	SRC/route/rr_graph_timing_params.c
	SRC/route/rr_graph_util.c
	SRC/route/segment_stats.c
	SRC/route/test.c
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

target_include_directories(Vpr PUBLIC 
	#${CMAKE_CURRENT_SOURCE_DIR}
	#printhandler/SRC/TIO_InputOutputHandlers 
	#libarchfpga/include 
	/usr/local/include
	/opt/local/include
	~/Downloads/tbb44_20150728oss/include
	SRC/util
	SRC/timing 
	SRC/pack 
	SRC/place 
	SRC/base 
	SRC/route 
	SRC/power 
)

target_compile_options(Vpr PRIVATE
	-Wno-write-strings -Wno-unused-variable -Wno-unused-private-field
	-DNO_GRAPHICS) 

target_link_libraries(Vpr
	ArchFpga
	PrintHandler
	Pcre
	)
