#include <iostream>
#include <fstream>
using std::cout; using std::cerr; using std::endl;
#include <vector>
using std::vector;
#include <stdexcept>
using std::domain_error;
#include <string>
using std::string;


#include <stdlib.h>

#include "hashmap.h"
#include "nautyinit.h"
#include "init.hpp"
#include "maingraph.hpp"
#include "output.hpp"

/*********************
 *  UTILITIES        *
 *********************/

// Writes the representation of a double value with 3 significant numbers to cout
void write_sci_repres(double number)
{
    int displ_exponent = 0;
	while (number >= 10) {
        number /= 10;
        ++displ_exponent;
    }
    std::streamsize prec = cout.precision();
	cout.precision(2);
    cout << number << "*10^" << displ_exponent;
	cout.precision(prec);
}

// Calculates the expected number of sampled subgraph
uint64 calc_expected_samples(double ntrees, const short & G_N, const double *prob)
{
    for (int i = 0; i != G_N; ++i)
    {
   	    ntrees *= prob[i];
	}
    uint64 SMPLS = (uint64)(ntrees+0.5);
	cout << " - Expected number of sampled subgraphs: ";
    write_sci_repres(ntrees);
    cout << endl;
    return SMPLS;
}

/*********************
 *  MAIN             *
 *********************/

int main(int argc, char **argv)
{
		
    // Welcome user
	cout << endl
	<< "------------------------------" << endl
	<< "          randESU 1.1         " << endl
	<< "           2004-05            " << endl
	<< "      Sebastian Wernicke      " << endl
	<< "        Florian Rasche        " << endl
	<< "------------------------------" << endl << endl;


    // Process flags and check correct usage
	cout << "1. Process Input Parameters and Initialize" << endl << endl;
	bool directed = true;
	bool fullenumeration = true;
    uint64 SMPLS = 1000000;
	uint64 TREESMPLS = 1000000; // Number of samples to estimate treesize
	short G_N = 4;
    char *inputfile = argv[1];
    string outputfile;
    double prob[] = { 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
    int num_r_nets = 0; // Number of random networks to generate
    short random_type = GLOBAL_CONST; // How to generate random network
    int num_exchanges = 3; // Number of exchanges per edge
    int num_tries = 3;     // Number of tries to find an exchange partner
    bool reest_size = false; // Estimate size for each generated graph
    long total_tries = 0, total_success = 0;
    string status_string, erase_string, empty_string;

    int ret = process_flags(argc, argv, SMPLS, TREESMPLS, G_N, 
		        directed, outputfile, prob, fullenumeration, num_r_nets, random_type, num_exchanges, num_tries, reest_size);
    if (ret != -1)
	  return ret;
	if (!directed) random_type = LOCAL_CONST; // any other type would be nonsense

    // Initialize Nauty
	cout << " - Initializing Nauty (c) Brendan McKay" << endl;
    const int G_M = (G_N + WORDSIZE - 1) / WORDSIZE;
    graph canon[MAXN * MAXM];
    graph g[MAXN * MAXM];
    for (int i = 0; i != G_N; ++i) {
		EMPTYSET(GRAPHROW(g, i, G_M), G_M);
    }
	int lab[MAXN], ptn[MAXN], orbits[MAXN];
    static DEFAULTOPTIONS(options);
    statsblk(stats);
    setword workspace[160 * MAXM];
    set *gv;
    options.writeautoms = FALSE;
    options.getcanon = TRUE;
    if (directed) {
		options.digraph = TRUE;
		options.invarproc = adjacencies;
		options.mininvarlevel = 1;
		options.maxinvarlevel = 10;
    }
    nauty_check(WORDSIZE, G_M, G_N, NAUTYVERSIONID);

	// Initialize Graph structure from file
    cout << endl << "2. Build Graph" << endl << endl;
    maingraph maing;
    try {
      maing = read_graph(inputfile, directed);
    } catch(domain_error e) {
		cerr << e.what() << endl;
		return 1;
    }


    // Build the arrays which store neighbour information
    alloc_graph_arrays(maing);
    build_graph(maing);
    // v_extension is used in the sampling and est_tree_size-functions
   	register long* v_extension = new long[maing.maxnumlargerneighbors*G_N];

    // Write general info into the outputfile
	std::ofstream outfile (outputfile.c_str());
    outfile << "Randesu 1.1 subgraph "
            << (fullenumeration ? "enumeration" : "sampling") << endl;
    outfile << "-------------------------------\n" << endl;
    outfile << "Network name: " << argv[1] << endl;
    outfile << "Network type: " << (maing.directed? "Directed" : "Undirected") << endl;
    outfile << "Number of nodes: " << maing.n << endl
            << "Number of edges: " << maing.m << endl;
    outfile << "Number of single edges: " << maing.num_dir_edges << endl
            << "Number of mutual edges: " << maing.num_undir_edges << endl << endl;
    outfile << "Algorithm: " << (fullenumeration ? "enumeration" : "sampling") << endl;
    if (!fullenumeration) {
        outfile << "Sampling parameters = {";
		for (int i = 0; i!= G_N; ++i)
			outfile << " " << prob[i];
		outfile << " }" << endl;
    }
    outfile << "Subgraph size: " << G_N << endl << endl;
    if (num_r_nets == 0)
    {
        outfile << "No random graphs were generated " << endl;
    } else {
    string type_str;
    switch (random_type) {
        case 0: type_str = "no respect to bidirectional edges"; break;
        case 1: type_str = "globally constant number of bidirectional edges";
                break;
        case 2: type_str = "locally constant number of bidirectional edges";
                break;
    }
    outfile << "Generated " << num_r_nets << " random networks" << endl
            << "   with " << type_str << endl
            << "   " << num_exchanges << " exchanges per edge and "
            <<  num_tries << " tries per edge" << endl << endl;
    }
    outfile.close(); // Close the file so it is not open during main algorithm.

    // Initializers for main algorithm 
    cout << endl << "3. Main Algorithm" << endl << endl;

	// Estimate number of subgraphs by random tree sampling
	uint64 APPROX_TREE_SIZE = est_tree_size(maing, v_extension, TREESMPLS, G_N);
	uint64 numtrees = 0x28F5C28F5C28F5CULL;

	if (TREESMPLS > 0) {
        double ntrees = (double(APPROX_TREE_SIZE) / double(TREESMPLS))*double(maing.n);
	    cout << " - Approximate number of subgraphs: ";
        write_sci_repres(ntrees);
        cout << endl
		 << "   (based on " << TREESMPLS << " samples)" << endl; 
    	numtrees = uint64(ntrees);
	    if (!fullenumeration) SMPLS = calc_expected_samples(ntrees, G_N, prob);
	}
    // For the Status bar: pass the number of expected subgraphs to sampling
    uint64 equiv100p = fullenumeration ? numtrees : SMPLS;

	// Start the sampling / enumeration of subgraphs
    cout << " - Executing subgraph "
	<< (fullenumeration ? "enumeration" : "sampling")
	<< " algorithm" << endl << endl;


	// Main enumeration / sampling loop 
    int total_num_nets = num_r_nets + 1; // We sample num_r_nets and the original graph
    // In this hashmap, the intermediate result is stored before entered into the result - hashmap.
    hash_map < graph64, uint64 > inter_result;
    // In this hashmap, the results of the sampling and the randomization are stored.
    hash_map < graph64, uint64* > result_graphs;
    uint64 *count_subgr = new uint64[total_num_nets];
    uint64 total_subgr = 0;
    double sampling_time = 0.0, random_time = 0.0;
    uint64 *current_array; // Abbrevation for the array currently updated

    // Sample the original graph and random graphs, if necessary.
    for (int nets_ctr = 0; nets_ctr < total_num_nets; ++nets_ctr){
       if (nets_ctr == 0){
           status_string = "Sampling original network: ";
       } else {
           status_string = "Sampling random network " + int_to_str(nets_ctr) +
                           " of " + int_to_str(num_r_nets) + ": ";
       }
       cout << status_string;
       cout.flush();
       // erase_string is as long as the status_string + the length of the statusbar
       erase_string.assign(status_string.size()+31, '\b');
       empty_string.assign(status_string.size()+31, ' ');

       // Reestimates tree-size if the user wishes.
       // and it is not the original network being sampled
       if (reest_size && nets_ctr != 0 && TREESMPLS > 0) {
           numtrees = est_tree_size(maing, v_extension, TREESMPLS, G_N) / TREESMPLS * maing.n;
	       if (!fullenumeration) SMPLS = calc_expected_samples(double(numtrees), G_N, prob);
           equiv100p = fullenumeration ? numtrees : SMPLS;
       }
       // Sample the graph
       // sampling is called with the Nauty-Vars, its own options,
       // true to show the status-bar and inter_result and count_subgr as result-parameters
       sampling_time += sampling(maing, v_extension, G_N, G_M,  fullenumeration, prob, equiv100p,
                                 true, canon, g, gv, lab, ptn, orbits, options,
                                 stats, workspace, inter_result, count_subgr[nets_ctr]);
       total_subgr += count_subgr[nets_ctr];
       // Write results into result_graphs-hashmap
       for (hash_map < graph64, uint64 >::const_iterator iter = inter_result.begin();
	        iter !=inter_result.end(); ++iter){
           if (result_graphs.find(iter->first) == result_graphs.end()){
               // create the array at the graph's hashmap position.
               current_array = (result_graphs[iter->first] = new uint64[total_num_nets]);
               // Initialize the array
               for (int i=0; i < total_num_nets; ++i)
                    current_array[i] = 0;
           } else {
               // set current_array as the array at the graphs map position
               current_array = result_graphs[iter->first];
           }
           // Write the number of subgraphs in the array
           current_array[nets_ctr] = iter->second;
       } // end for
       inter_result.clear();       // Empty the intermediate result hashmap.

       // Randomize the graph
       if (num_r_nets != 0)
         random_time += randomize_graph(maing, random_type, num_exchanges, num_tries, total_tries, total_success);
       // Erase the Statusbar line:
       cout << erase_string << empty_string << erase_string;
    } // end for

    // Free memory
    maing.edges.clear();
    maing.num_neighbours.clear();
    maing.num_undir.clear();
    for (vertex i = 0; i != maing.n; ++i) {
      delete[] maing.neighbours[i] ;
    }
    delete[] maing.neighbours;
    delete[] maing.num_larger_neighbours;
    delete[] maing.v_util;
    delete[] v_extension;

	// Output execution statistics
    cout << "4. Results" << endl << endl;
    std::streamsize prec = cout.precision();
    cout.precision(2);
    cout << " - Overall " << (fullenumeration ? "enumeration" : "sampling") 
		 << " time: " << sampling_time << " seconds" << endl
	     << "   (equals "
	     << 1000000.0 * sampling_time / (double) total_subgr;
	if (1000000.0 * 0.01 / (double) total_subgr > 0.1 ) {
	     cout << " +/- " 
	          << 1000000.0 * 0.01 / (double) total_subgr ;
	}
	cout << " microseconds per subgraph)" << endl;
    if (num_r_nets != 0){
        cout << " - Overall randomization time: " << random_time << " seconds"<< endl;
    }
    cout << " - Number of "
	     << (fullenumeration ? "enumerated" : "sampled") 
	     << " subgraphs in original network: " << count_subgr[0] << endl;
    if(num_r_nets != 0){
        cout << " - Number of "
	         << (fullenumeration ? "enumerated" : "sampled")
	         << " subgraphs in random networks: "
             << (total_subgr-count_subgr[0]) << endl;
        cout << " - Number of "
	         << (fullenumeration ? "enumerated" : "sampled")
	         << " subgraphs in all networks: " << total_subgr << endl;
    }
    cout << " - Number of subgraph classes in all networks: " << result_graphs.size() << endl;
    cout.precision(prec);

	// Output results in file
    outfile.open(outputfile.c_str(), std::ofstream::out | std::ofstream::app);
    outfile << count_subgr[0] << " subgraphs were "
            << (fullenumeration ? "enumerated" : "sampled")
            << " in the original network." << endl;
    if (num_r_nets != 0) {
        outfile << (total_subgr-count_subgr[0]) << " subgraphs were "
                << (fullenumeration ? "enumerated" : "sampled")
                << " in the random networks." << endl;
        outfile << total_subgr << " subgraphs were "
                << (fullenumeration ? "enumerated" : "sampled")
                << " in all networks." << endl << endl;

        outfile << "For the random networks: " << total_tries << " tries were made, "
                << total_success << " were successful." << endl;
        outfile << "Randomization took " << random_time << " seconds." << endl;
    } else {
        outfile << endl;
    }
    outfile << (fullenumeration ? "Enumeration took " : "Sampling took ")
            << sampling_time << " seconds.\n\n" << endl;
    // Output the results table
	pretty_output(result_graphs, G_N, count_subgr, num_r_nets, outfile);
    outfile.close();

    // Free Memory of result storage
    for (hash_map < graph64, uint64* >::const_iterator iter = result_graphs.begin();
         iter != result_graphs.end(); ++iter){
        delete[] iter->second;
    }
    result_graphs.clear();
    delete[] count_subgr;

	// Done.
	cout << " - Results have been written to \'" 
		 << outputfile << "\'" << endl;
	return 0;
}
