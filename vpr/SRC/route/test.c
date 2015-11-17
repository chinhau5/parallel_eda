#include <stdio.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/numeric/interval.hpp>
#include <boost/timer/timer.hpp>
#include <vector>
extern "C" {
#include <igraph/igraph.h>
}
#include <stdlib.h>
#include "vpr_types.h"
#include "route.h"
#include "graph.h"


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<int, 2, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, int> value;

void test_rtree()
{
	std::vector<value> bulk;
	box b1(point(0, 0), point(2, 2));
	/*box b(point(), point());*/
	/*box b(point(), point());*/
	// insert new value
	bulk.push_back(std::make_pair(b1, bulk.size()));
	/*bulk.push_back(std::make_pair(b2, bulk.size()));*/
	/*bulk.push_back(std::make_pair(b3, bulk.size()));*/

	bgi::rtree< value, bgi::rstar<16> > t(bulk);

	box query_box(point(2,2), point(3, 3));
	std::vector<value> result_s;
	t.query(bgi::intersects(query_box), std::back_inserter(result_s));

	for (const auto &v : result_s) {
		printf("overlapping: %d\n", v.second);
	}
}

void print_vector(igraph_vector_t *v) {
  long int i, n=igraph_vector_size(v);
  for (i=0; i<n; i++) {
    printf(" %li", (long int) VECTOR(*v)[i]);
  }
  printf("\n");
}

void warning_handler_ignore(const char* reason,const char* file,int line,int e) {
}

int schedule_nets_maximum_independent_set(const RRGraph &rr_g)
{  
  igraph_t g;
  igraph_vector_ptr_t result;
  long int i, j, n;
  igraph_integer_t alpha;
  const int params[] = {4, -1, 2, 2, 0, 0, -1, -1};
  igraph_vector_t edges;
 
  igraph_set_warning_handler(warning_handler_ignore);
  igraph_vector_init(&edges, num_edges(rr_g)*2);
  igraph_vector_ptr_init(&result, 0);

  int edge = 0;
  for_all_edges(rr_g, [&rr_g, &edges, &edge] (const RREdge &e) -> void {
		  int from = id(get_source(rr_g, e));
		  int to = id(get_source(rr_g, e));
		  VECTOR(edges)[edge] = from;
		  VECTOR(edges)[edge+1] = to;
		  edge += 2;
		  });

  igraph_create(&g, &edges, num_vertices(rr_g), 1);
  igraph_largest_independent_vertex_sets(&g, &result);
  n = igraph_vector_ptr_size(&result);
  printf("%ld largest independent sets found\n", (long)n);
  for (i=0; i<n; i++) {
	igraph_vector_t* v = (igraph_vector_t *)igraph_vector_ptr_e(&result,i);
	int set_size = igraph_vector_size(v);
	printf("Independent set %d size: %d\n", i, set_size);
	print_vector((igraph_vector_t*)v);
	igraph_vector_destroy(v);
	free(v);
  }
  igraph_vector_ptr_destroy(&result);

  /*igraph_tree(&g, 5, 2, IGRAPH_TREE_OUT);*/
  /*for (j=0; j<sizeof(params)/(2*sizeof(params[0])); j++) {*/
    /*if (params[2*j+1] != 0) {*/
      /*igraph_independent_vertex_sets(&g, &result, params[2*j], params[2*j+1]);*/
    /*} else {*/
      /*igraph_largest_independent_vertex_sets(&g, &result);*/
    /*}*/
    /*n = igraph_vector_ptr_size(&result);*/
    /*printf("%ld independent sets found\n", (long)n);*/
    /*for (i=0; i<n; i++) {*/
      /*igraph_vector_t* v;*/
      /*v=(igraigraph_vector_ptr_e(&result,i);*/
      /*print_vector((igraph_vector_t*)v);*/
      /*igraph_vector_destroy(v);*/
      /*free(v);*/
    /*}*/
  /*}*/
  /*igraph_destroy(&g);*/

  /*igraph_tree(&g, 10, 2, IGRAPH_TREE_OUT);*/
  /*igraph_maximal_independent_vertex_sets(&g, &result);*/
  /*n = igraph_vector_ptr_size(&result);*/
  /*printf("%ld maximal independent sets found\n", (long)n);*/
  /*for (i=0; i<n; i++) {*/
    /*igraph_vector_t* v;*/
    /*v=igraph_vector_ptr_e(&result,i);*/
    /*print_vector((igraph_vector_t*)v);*/
    /*igraph_vector_destroy(v);*/
    /*free(v);*/
  /*}*/
  /*igraph_vector_ptr_destroy(&result);*/

  /*igraph_independence_number(&g, &alpha);*/
  /*printf("alpha=%ld\n", (long)alpha);*/

  igraph_destroy(&g);

  return 0;
}
