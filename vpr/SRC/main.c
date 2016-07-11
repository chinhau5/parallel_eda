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

#include "route.h"

#include <sched.h>
#include <numaif.h>
#include <numa.h>
#include <unistd.h>

/*#include "vt_user.h"*/

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

void print_mem_bind()
{
	int num_possible_nodes = numa_num_possible_nodes();

	assert(num_possible_nodes % sizeof(unsigned long) == 0);

	int nodemask_size = num_possible_nodes / sizeof(unsigned long);

	printf("nodemask size %d\n", nodemask_size);

	int mode;
	unsigned long *nodemask = (unsigned long *)malloc(nodemask_size);
	if (get_mempolicy(&mode, nodemask, num_possible_nodes, 0, 0) != 0) {
		perror("get_mempolicy error");
		exit(-1);
	}

	printf("mempolicy mode: ");
	switch (mode) {
		case MPOL_DEFAULT:
			printf("default\n");
			break;
		case MPOL_BIND:
			printf("bind\n");
			break;
		case MPOL_INTERLEAVE:
			printf("interleave\n");
			break;
		case MPOL_PREFERRED:
			printf("preferred\n");
			break;
		default:
			printf("unknown\n");
			break;
	}

	printf("mem_bind: ");
	int curnode = 0;
	for (int i = 0; i < 1 && curnode < num_possible_nodes; ++i) {
		unsigned long val = nodemask[i];
		for (int j = 0; j < sizeof(unsigned long)*8 && curnode < num_possible_nodes; ++j) {	
			if (val & 1) {
				printf("%d ", curnode);
			}
			val >>= 1;
			++curnode;
		}
	}
	printf("\n");
}

void print_mems_allowed()
{
	struct bitmask *bm = numa_get_mems_allowed();
	printf("mems_allowed: \n");
	int i;
	for (i = 0; i < numa_num_possible_nodes(); ++i) {
		if (numa_bitmask_isbitset(bm, i)) {
			printf("%d ");
		}
	}
	printf("\n");
}

void get_sched_bind(vector<int> &cpuset, set<int> &nodeset)
{
	int num_cpus = numa_num_configured_cpus();
	printf("num cpus: %d\n", num_cpus);
	cpu_set_t *cpumask = CPU_ALLOC(num_cpus);
	printf("cpu_set_t alloc size: %d\n", CPU_ALLOC_SIZE(num_cpus));
	if (sched_getaffinity(0, CPU_ALLOC_SIZE(num_cpus), cpumask)) {
		perror("Error getting affinity\n");
		exit(-1);
	}

	for (int i = 0; i < num_cpus; ++i) {
		if (CPU_ISSET(i, cpumask)) {
			cpuset.push_back(i);
			nodeset.insert(numa_node_of_cpu(i));
		}
	}
}

void print_sched_bind(const vector<int> &cpuset, const set<int> &nodeset)
{
	printf("cpu affinity: ");
	for (const auto &cpu : cpuset) {
		printf("%d ", cpu);
	}
	printf("\n");

	printf("node affinity: ");
	for (const auto &n : nodeset) {
		printf("%d ", n);
	}
	printf("\n");
}

int main(int argc, char **argv) {
	t_options Options;
	t_arch Arch;
	t_vpr_setup vpr_setup;
	clock_t entire_flow_begin,entire_flow_end;

	vector<int> cpuset;
	set<int> nodeset;

	get_sched_bind(cpuset, nodeset);

	print_sched_bind(cpuset, nodeset);
	print_mem_bind();

	bitmask *bind_to = numa_allocate_nodemask();
	numa_bitmask_clearall(bind_to);
	for (const auto &n: nodeset) {
		numa_bitmask_setbit(bind_to, n);
	}
	numa_bind(bind_to);

	bitmask *from = numa_allocate_nodemask();
	numa_bitmask_setall(from);

	numa_migrate_pages(0, from, bind_to);

	/*int net_disp[] = {*/
		/*offsetof(net_t, vpr_id),*/

		/*[> source <]*/
		/*offsetof(net_t, source)+offsetof(source_t, rr_node),*/
		/*offsetof(net_t, source)+offsetof(source_t, x),*/
		/*offsetof(net_t, source)+offsetof(source_t, y),*/
		
		/*[> bounding box <]*/
		/*offsetof(net_t, bounding_box)+offsetof(bounding_box_t, xmin),*/
		/*offsetof(net_t, bounding_box)+offsetof(bounding_box_t, xmax),*/
		/*offsetof(net_t, bounding_box)+offsetof(bounding_box_t, ymin),*/
		/*offsetof(net_t, bounding_box)+offsetof(bounding_box_t, ymax)*/
	/*};*/

	/*for (int i = 0; i < 7; ++i) {*/
		/*printf("Offset is %X\n", net_disp[i]);*/
	/*}*/
	/*return 0;*/

	/*MPI_Init(&argc, &argv);*/

	/*mtrace();*/

	/*VT_USER_START("test");*/

	program_start = std::chrono::high_resolution_clock::now(); 

	entire_flow_begin = clock();

	/*VT_USER_END("test");*/

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




