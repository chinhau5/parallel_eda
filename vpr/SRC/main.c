/**
 VPR is a CAD tool used to conduct FPGA architecture exploration.  It takes, as input, a technology-mapped netlist and a description of the FPGA architecture being investigated.  
 VPR then generates a packed, placed, and routed FPGA (in .net, .place, and .route files respectively) that implements the input netlist.
 
 This file is where VPR starts execution.

 Key files in VPR:
 1.  libarchfpga/physical_types.h - Data structures that define the properties of the FPGA architecture
 2.  vpr_types.h - Very major file that defines the core data structures used in VPR.  This includes detailed architecture information, user netlist data structures, and data structures that describe the mapping between those two.
 3.  globals.h - Defines the global variables used by VPR.
 */


#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
/*#include <glog/logging.h>*/
#include <zlog.h>
#include <chrono>
#include <mpi.h>
/*#include <mcheck.h>*/
/*#include "tbb/task_scheduler_init.h"*/
#include "vpr_api.h"

/**
 * VPR program
 * Generate FPGA architecture given architecture description
 * Pack, place, and route circuit into FPGA architecture
 * Electrical timing analysis on results
 *
 * Overall steps
 * 1.  Initialization
 * 2.  Pack
 * 3.  Place-and-route and timing analysis
 * 4.  Clean up
 */

void test_dfs();
void init_parallel_route_logging();
void init_advanced_parallel_route_logging();

std::chrono::time_point<std::chrono::high_resolution_clock> program_start;

char *s_circuit_name = nullptr;

int main(int argc, char **argv) {
	t_options Options;
	t_arch Arch;
	t_vpr_setup vpr_setup;
	clock_t entire_flow_begin,entire_flow_end;

	MPI_Init(&argc, &argv);

	/*mtrace();*/

	program_start = std::chrono::high_resolution_clock::now(); 

	entire_flow_begin = clock();

	if (dzlog_init("log.conf", "default") == -1) {
		printf("failed to init zlog\n");
		return -1;
	}

	/*init_parallel_route_logging();*/
	/*init_advanced_parallel_route_logging();*/

	/*test_dfs();*/
	/*return 0;*/

	/*google::InitGoogleLogging(argv[0]);*/
	/*LOG(INFO) << "test";*/

	/* Read options, architecture, and circuit netlist */
	vpr_init(argc, argv, &Options, &vpr_setup, &Arch);

	if (setvbuf(stdout, NULL, _IONBF, 0)) {
		exit(-1);
	}

	/* If the user requests packing, do packing */
	if (vpr_setup.PackerOpts.doPacking) {
		vpr_pack(vpr_setup, Arch);
	}

	if (vpr_setup.PlacerOpts.doPlacement || vpr_setup.RouterOpts.doRouting) {
		vpr_init_pre_place_and_route(vpr_setup, Arch);
		vpr_place_and_route(vpr_setup, Arch);
#if 0
		if(vpr_setup.RouterOpts.doRouting) {
			vpr_resync_post_route_netlist_to_TI_CLAY_v1_architecture(&Arch);
		}
#endif
	}

	if (vpr_setup.PowerOpts.do_power) {
		vpr_power_estimation(vpr_setup, Arch);
	}
	
	entire_flow_end = clock();
	
	#ifdef CLOCKS_PER_SEC
		vpr_printf(TIO_MESSAGE_INFO, "The entire flow of VPR took %g seconds.\n", (float)(entire_flow_end - entire_flow_begin) / CLOCKS_PER_SEC);
	#else
		vpr_printf(TIO_MESSAGE_INFO, "The entire flow of VPR took %g seconds.\n", (float)(entire_flow_end - entire_flow_begin) / CLK_PER_SEC);
	#endif
	
	/* free data structures */
	vpr_free_all(Arch, Options, vpr_setup);

	/* Return 0 to single success to scripts */
	return 0;
}




