#ifndef MAINGRAPH_HPP
#define MAINGRAPH_HPP

#include <iostream>
using std::cout; using std::cerr; using std::endl;
#include <stdexcept>
using std::domain_error;
#include <vector>
using std::vector;
#include <string>
using std::string;

extern "C" {
#include <cstdio>
#include <ctime>
}

#include "random.hpp"
#include "hashmap.h"
#include "nautyinit.h"
#include "graph64.hpp"


typedef struct
{
    hash_map < edge, edgetype > edges;
    vector < long > num_neighbours;
    vector < long > num_undir; // stores the number of undirected edges from a vertex
    long n, m; // Number vertices, Number edges
    long num_nodes;
    long num_lonely_nodes;
    long num_dir_edges;
    long num_undir_edges;
    bool directed;

    vertex **neighbours;
	long *num_larger_neighbours;
	long maxnumlargerneighbors;
    unsigned long *v_util;

}  maingraph;

// These are constants for the random types
const short NO_REGARD = 0;    // No regard to bidirectional edges
const short GLOBAL_CONST = 1; // Keeping bidir edges globally const
const short LOCAL_CONST = 2;  // Keeping bidir edges const for each vertex

// Status bar functions:
void write_statusbar(int & index, const uint64 & equiv100p, const uint64 & status);
void finish_statusbar(int & index);

// Maingraph functions:
maingraph read_graph(char *filename, bool directed);
void alloc_graph_arrays(maingraph & maing);
void build_graph(maingraph & maing);
uint64 est_tree_size(const maingraph & g, long* v_extension, uint64 TREESMPLS, short G_N);
double sampling(const maingraph & maing, long* v_extension, short G_N, const int G_M,  bool fullenumeration, const double* prob,
                const uint64 equiv100p, const bool write_status, graph* canon, graph* g, set* gv, int* lab, int* ptn, int* orbits, optionblk & options,
                  statsblk & stats, setword* workspace, hash_map < graph64, uint64 > & result_graphs, uint64 & count_subgr);
double randomize_graph(maingraph & maing, short random_type,
                     int num_exchanges, int num_tries, long & total_tries,
                     long & total_success);

#endif
