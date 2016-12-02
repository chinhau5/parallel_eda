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
#include "config.h"

#ifdef __linux__
#include <sched.h>
#ifdef WITH_NUMA
#include <numaif.h>
#include <numa.h>
#endif
#include <unistd.h>
#include <malloc.h>
#include <sys/resource.h>
#endif

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

#if defined(__linux__) && defined(WITH_NUMA)

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

void print_mem_map()
{
	int pid = getpid();
	char buffer[512];
	sprintf(buffer,"/proc/%d/maps", pid);
	FILE *fd = fopen(buffer, "r");
	while (!feof(fd)) {

		fgets(buffer, 512, fd);
		printf(buffer);
	}
	fclose(fd);
}

void print_env()
{
	extern char **environ;
		
	printf("Environment start\n");
	for (char **cur = environ; *cur != 0; ++cur) {
		printf("%s\n", *cur);
	}
	printf("Environment end\n");
}

typedef struct mem_map_entry_t {
	unsigned long start;
	unsigned long end;
	unsigned long file_offset;
	char image_name[256];
} mem_map_entry_t;

void parse_mem_map(int pid, vector<mem_map_entry_t> &entries)
{
	char buffer[512];

	sprintf(buffer,"/proc/%d/maps", pid);

	FILE *fd = fopen(buffer, "r");

	while (!feof(fd)) {
		fgets(buffer, 512, fd);

		mem_map_entry_t entry;
		sscanf(buffer, "%lX-%lX %*s %lX %*s %*s %s", &entry.start, &entry.end, &entry.file_offset, entry.image_name);
		/*printf("%lX-%lX %lX %s\n", entry.start, entry.end, entry.file_offset, entry.image_name);*/

		entries.push_back(entry);
	}

	fclose(fd);
}

const mem_map_entry_t *find_sym(const vector<mem_map_entry_t> &entries, unsigned long addr)
{
	for (const auto &entry : entries) {
		if (addr >= entry.start && addr <= entry.end) {
			return &entry;
		}
	}
	return nullptr;
}

string find_func_name(const mem_map_entry_t *entry, unsigned long addr)
{
	assert(entry->file_offset == 0);
	unsigned long offset = addr-entry->start;
	char cmd[256];
	sprintf(cmd, "addr2line -e %s -f 0x%lX", entry->image_name, offset);
	FILE *pipe = popen(cmd, "r");
	int num_lines = 0;
	string func_name = "";
	while (!feof(pipe)) {
		char buffer[256];
		fgets(buffer, 256, pipe);
		printf("addr2line %d: %s\n", num_lines, buffer);
		if (num_lines == 0) {
			func_name = buffer;	
		}
		++num_lines;
	}
	pclose(pipe);

	return func_name;
}

#endif

/*extern void *(*__malloc_hook)(size_t size, const void *caller);*/
void print_context(int pid, int rank)
{
#if defined(__linux__)
	/*pid_t pid = getpid();*/

	/*char buffer[256];*/
	/*sprintf(buffer, "/proc/%d/status", pid);*/

	/*FILE *stat = fopen(buffer, "r");*/
	/*while (!feof(stat)) {*/
		/*int num_read = fread(buffer, 1, 256, stat);*/
		/*buffer[num_read] = 0;*/
		/*printf("%s", buffer);*/
	/*}*/

	printf("[%d,%d] LSB_BIND_CPU_LIST: %s\n", pid, rank, getenv("LSB_BIND_CPU_LIST"));

	cpu_set_t cpuset;
	assert(sched_getaffinity(0, sizeof(cpu_set_t), &cpuset) == 0);

	printf("[%d,%d] cpu affinity [count %d max %d]: ", pid, rank, CPU_COUNT(&cpuset), CPU_SETSIZE);
	for (int i = 0; i < CPU_SETSIZE; ++i) {
		if (CPU_ISSET(i, &cpuset)) {
			printf("%d ", i);
		}
	}
	printf("\n");

	char hostname[256];
	gethostname(hostname, 256);
	printf("[%d,%d] hostname: %s\n", pid, rank, hostname);
#endif
}

int main(int argc, char **argv) {
	t_options Options;
	t_arch Arch;
	t_vpr_setup vpr_setup;
	clock_t entire_flow_begin,entire_flow_end;
	
	print_context(getpid(), -1);

#if defined(__linux__) && defined(WITH_NUMA)
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

	vector<mem_map_entry_t> entries_before_init;
	parse_mem_map(getpid(), entries_before_init);
	unsigned long malloc_hook_addr = (unsigned long)__malloc_hook;
	const mem_map_entry_t *entry = find_sym(entries_before_init, malloc_hook_addr);
	if (entry) {
		string func_name = find_func_name(entry, malloc_hook_addr);
		printf("malloc_hook before mpi_init: %s (%lX)\n", func_name.c_str(), __malloc_hook);
	} else {
		printf("malloc_hook before mpi_init not found: %lX\n", __malloc_hook);
	}
#endif

#if defined(__linux__) 
	struct rlimit lim;
	lim.rlim_cur = 10240*1024;
	lim.rlim_max = 10240*1024;

	if (setrlimit(RLIMIT_STACK, &lim) != 0) {
		printf("failed to set limit\n");
		return -1;
	}
#endif

#if defined(__linux__) 
	cpu_set_t cpuset;
	CPU_ZERO(&cpuset);
	CPU_SET(0, &cpuset);

	/*assert(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) == 0);*/
#endif

#ifdef VPR_MPI
	MPI_Init(&argc, &argv);
#endif

#if defined(__linux__) && defined(WITH_NUMA)
	/*print_mem_map();*/
	/*print_env();*/

	vector<mem_map_entry_t> entries_after_init;
	parse_mem_map(getpid(), entries_after_init);
	malloc_hook_addr = (unsigned long)__malloc_hook;
	entry = find_sym(entries_after_init, malloc_hook_addr);
	if (entry) {
		string func_name = find_func_name(entry, malloc_hook_addr);
		printf("malloc_hook after mpi_init: %s (%lX)\n", func_name.c_str(), __malloc_hook);
	} else {
		printf("malloc_hook after mpi_init not found: %lX\n", __malloc_hook);
	}
#endif

	int rank = 0;
#ifdef VPR_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

	print_context(getpid(), rank);

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




