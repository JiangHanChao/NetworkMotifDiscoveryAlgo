#include <string>
using std::string;
#include <sstream>
using std::stringstream;
#include <fstream>

#include "graph64.hpp"

int process_flags(int argc, char **argv, uint64 & SMPLS, uint64 & TREESMPLS,
	      short &G_N, bool & directed, string & outputfile,
	      double *prob, bool & fullenumeration, int & num_r_nets, short & random_type,
          int & num_exchanges, int & num_tries, bool & reest_size);
