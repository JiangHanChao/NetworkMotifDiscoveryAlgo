#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "graph64.hpp"
#include "maingraph.hpp"
#include "hashmap.h"

struct subgr_res {
    graph64 id;
    double freq;
    double rand_mean;
    double rand_sd;
    double p_value;
    double z_score;
};

inline bool compare(const subgr_res & a, const subgr_res & b)
{
    return a.freq > b.freq;
}

// Converts an int into string
std::string int_to_str(const int & i);

void
pretty_output(hash_map < graph64, uint64* > & res_graphs, short G_N,
              uint64* count_subgr, int num_r_nets, std::ofstream & outfile);

#endif
