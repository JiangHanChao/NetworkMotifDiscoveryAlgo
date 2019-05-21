#include "init.hpp"

int
process_flags(int argc, char **argv, uint64 & SMPLS, uint64 & TREESMPLS,
	      short &G_N, bool & directed, string & outputfile,
	      double *prob, bool & fullenumeration, int & num_r_nets, short & random_type,
          int & num_exchanges, int & num_tries, bool & reest_size)
{
    if (argc < 2) {
	cerr << "ERROR: Usage is \'" << argv[0] << " <infile> [options]\'"
	    << endl;
	cerr << "       Call \'" << argv[0] << " --help\' for more information." << endl;
	return 1;
    }
    // Initialize
    bool subgraphsize_set = false,
	outfilename_set = false,
	directed_set = false,
	numsamples_set = false,
	trav_set = false,
	tree_set = false,
    random_set = false,
    random_t_set = false,
    random_ex_set = false,
    random_tr_set = false,
    reest_size_set = false;

    double p;
    directed = true;
    outputfile = string(argv[1]) + ".OUT";

    for (int i = 1; i != argc; ++i) {

	string s = string(argv[i]);
	if (s == "-s") {
	    ++i;
	    if (argc == i) {
		cerr << "ERROR: Flag \'-s\' expects number to follow." <<
		    endl;
		return 1;
	    }
	    string t = string(argv[i]);
	    stringstream ss(t);
	    ss >> G_N;
	    if (ss.fail() || G_N < 1 || G_N > 8) {
		cerr <<
		    "ERROR: Expected number between 2 and 8 to follow flag \'-s\'"
		    << endl;
		return 1;
	    }
	    subgraphsize_set = true;
	} else if (s == "-t") {
	    ++i;
	    if (argc == i) {
		cerr << "ERROR: Flag \'-t\' expects number to follow." <<
		    endl;
		return 1;
	    }
	    string t = string(argv[i]);
	    stringstream ss(t);
	    ss >> TREESMPLS;
	    if (ss.fail() || TREESMPLS < 0 || TREESMPLS > 10000000) {
		cerr <<
		    "ERROR: Expected number between 0 and 1e7 to follow flag \'-t\'"
		    << endl;
		return 1;
	    }
	    tree_set = true;

	} else if (s == "-o" || s == "-f") {
	    ++i;
	    if (argc == i) {
		cerr << "ERROR: Flag \'" << s <<
		    "\' expects filename to follow." << endl;
		return 1;
	    }
	    string t = string(argv[i]);
	    stringstream ss(t);
	    ss >> outputfile;
	    if (ss.fail()) {
		cerr << "ERROR: Expected filename to follow \'" << s <<
		    "\'" << endl;
		return 1;
	    }
	    outfilename_set = true;
	} else if (s == "-prob") {  //TODO
		fullenumeration = false;
	    ++i;
		if (argc == i) {
			cerr << "ERROR: Flag \'" << s <<
				"\' expects numbers between 0.0 and 1.0 to follow." <<
				endl;
			return 1;
		}
		int j = i;
		do
		{
			string t = string(argv[j]);
			stringstream ss(t);
			ss >> p;
			if (p < 0 || p > 1) {
				cerr << "ERROR: Expected numbers between 0.0 and 1.0 to follow"
					 << "\'" << s << "\'" << endl;
				return 1;
			}		
			if (ss.fail()) {
				if (j==i) {
					cerr << "ERROR: Expected numbers between 0.0 and 1.0 to follow"
						 << "\'" << s << "\'" << endl;
					return 1;
				} else {
					break;
				}
			}	
			if (j-i < 8)
			{  prob[j-i] = p;	}
			++j;
		} while (j!=argc);
		i = j-1;
	    trav_set = true;
	} else if (s == "-h" || s == "--help" || s == "-help") {
	    cout << "Usage of the program is \'" << argv[0]
		     << " <infile> [options]\'" << endl << endl;;
	    cout << "Recognized option flags are:" << endl;
	    cout <<	"  -s <int>        size of subgraphs to sample (2-8)    " <<	endl;
	    cout <<	"  -o <string>     name of the output file              " <<	endl;
	    cout <<	"  -f <string>     equivalent to -o                     " <<	endl;
	    cout <<	"  -directed       input graph is directed              " <<	endl;
	    cout <<	"  -d              equivalent to -directed              " <<	endl;
	    cout <<	"  -undirected     input graph is undirected            " <<	endl;
	    cout <<	"  -nondirected    equivalent to -undirected            " <<	endl;
	    cout <<	"  -nd             equivalent to -undirected            " <<	endl;
	    cout <<	"  -ud             equivalent to -undirected            " <<	endl;
	    cout <<	"  -t <int>        number of samples to estimate total  " <<    endl;
		cout << "                  number of subgraphs                  " <<    endl;
	    cout <<	"  -prob <doubles> activates sampling mode and sets treenode" << endl;
		cout << "                  acceptance probabilities (0..1)" <<	endl;
	    cout <<	"  -fullenum       full subgraph enumeration mode       " <<	endl;
	    cout <<	"  -r              Number of random networks to generate" <<	endl;
	    cout <<	"  -rt             Type of random subgraphs (0,1,2)     " <<	endl;
	    cout <<	"                  0 = Regardless of the number of bidir. edges." <<	endl;
	    cout <<	"                  1 = Keeping the bidir. edges globally constant" <<	endl;
	    cout <<	"                  2 = Keeping the bidir. edges constant for each node" <<	endl;
	    cout <<	"  -rn             Number of exchanges per edge         " <<	endl;
	    cout <<	"  -rs             Number of tries to search exchange partner" <<	endl;
	    cout <<	"  -re             Reestimate the number of subgraphs   " <<    endl;
        cout << "                  for each random graph                " <<	endl;
	    cout <<	"  -h              display this help                    " <<	endl;
	    cout <<	"  -help           equivalent to -h                     " <<	endl;
	    cout << "  --help          equivalent to -h             " <<	endl << endl;
	    cout << "Unspecified options imply default values.     " <<	endl << endl;
	    return 0;
	} else if (s == "-directed" || s == "-d") {
	    directed = true;
	    directed_set = true;
	} else if (s == "-fullenum" || s == "-fullenumeration") {
	    fullenumeration = true;
	    trav_set = true;
	} else if (s == "-undirected" || s == "-nondirected"
		   || s == "-nd" || s == "-ud") {
	    directed = false;
	    directed_set = true;
	} else if (s == "-r") {
	    ++i;
	    if (argc == i) {
		cerr << "ERROR: Flag \'-r\' expects number to follow." <<
		    endl;
		return 1;
	    }
	    string t = string(argv[i]);
	    stringstream ss(t);
	    ss >> num_r_nets;
	    if (ss.fail() || num_r_nets < 0 || num_r_nets > 10000000) {
		cerr <<
		    "ERROR: Expected number between 0 and 1e7 to follow flag \'-r\'"
		    << endl;
		return 1;
	    }
	    random_set = true;
	} else if (s == "-rt") {
	    ++i;
	    if (argc == i) {
		cerr << "ERROR: Flag \'-rt\' expects number to follow." <<
		    endl;
		return 1;
	    }
	    string t = string(argv[i]);
	    stringstream ss(t);
	    ss >> random_type;
	    if (ss.fail() || random_type < 0 || random_type > 2) {
		cerr <<
		    "ERROR: Expected number between 0 and 2 to follow flag \'-rt\'"
		    << endl;
		return 1;
	    }
	    random_t_set = true;
	} else if (s == "-rn") {
	    ++i;
	    if (argc == i) {
		cerr << "ERROR: Flag \'-rn\' expects number to follow." <<
		    endl;
		return 1;
	    }
	    string t = string(argv[i]);
	    stringstream ss(t);
	    ss >> num_exchanges;
	    if (ss.fail() || num_exchanges <= 0) {
		cerr <<
		    "ERROR: Expected number bigger than 0 to follow flag \'-rn\'"
		    << endl;
		return 1;
	    }
	    random_ex_set = true;
	} else if (s == "-rs") {
	    ++i;
	    if (argc == i) {
		cerr << "ERROR: Flag \'-rs\' expects number to follow." <<
		    endl;
		return 1;
	    }
	    string t = string(argv[i]);
	    stringstream ss(t);
	    ss >> num_tries;
	    if (ss.fail() || num_tries <= 0) {
		cerr <<
		    "ERROR: Expected number bigger than 0 to follow flag \'-rs\'"
		    << endl;
		return 1;
	    }
	    random_tr_set = true;
	} else if (s == "-re") {
	    reest_size = true;
	    reest_size_set = true;
	} else if (i != 1) {
	    cerr << "ERROR: Unrecognized flag \'" << s << "\'." << endl
		<< "       Call \'" << argv[0] << " --help\' for more "
		<< "information." << endl;
	    return 1;
	}
    } // end for

    cout << " - Algorithm: "
	     << (fullenumeration ? "enumeration" : "sampling") 
         <<	((!trav_set) ? " (default)" : "") << endl;
    if (!fullenumeration) {
		cout << "   sampling parameters = {";
		std::streamsize prec = cout.precision();
		cout.precision(2);
		for (int i = 0; i!= G_N; ++i)
			{ cout << " " << prob[i]; }
		cout.precision(prec);
		cout << " }" << endl;
	}
    cout << " - Subgraph size = " << G_N << ((!subgraphsize_set) ?
    				     " (default value)" : "") << endl;
	if (TREESMPLS == 0)
	{
	cout << " - Estimation of total subgraph number turned off" << endl
		 << "   [progress indicator will not work]" << endl;
	} else {
	cout << " - Samples to estimate no. of subgraphs = " << TREESMPLS
	     << ((!tree_set) ? " (default value)" : "") << endl;
	}
    cout << " - Input graph is " << (directed ? "" : "un")
		 << "directed " << ((!directed_set) ? "(default)" : "") << endl;
    if (num_r_nets == 0)
    {
    cout << " - No random graphs are generated "
         << ((!random_set) ? "(default)" : "") << endl;
    } else {
    string type_str;
    switch (random_type) {
        case 0: type_str = "no respect to bidirectional edges"; break;
        case 1: type_str = "globally constant number of bidirectional edges";
                break;
        case 2: type_str = "locally constant number of bidirectional edges";
                break;
    }
    cout << " - Generate " << num_r_nets << " random networks" << endl
         << "   with " << type_str << ((!random_t_set) ? " (default)," : ",") <<endl
         << "   " << num_exchanges << " exchanges per edge"
         << ((!random_ex_set) ? " (default), " : ", ") << num_tries << " tries per edge"
         << ((!random_tr_set) ? " (default)" : "") << endl
         << (reest_size ? "   and reestimation of subgraph number\n" : "") << endl;
    }
    cout << " - Name of output file = \'" << outputfile << "\'"
	     <<	((!outfilename_set) ? " (default value)" : "") << endl;
    return -1;
}
