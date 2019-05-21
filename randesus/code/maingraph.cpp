#include "maingraph.hpp"

/*********************
 *  UTILITIES        *
 *********************/

// Writes the Status bar according to actual status
void write_statusbar(int & index, const uint64 & equiv100p, const uint64 & status)
{
    while (status > (equiv100p * index * 5 / 100) && index < 20) {
        if (index % 5 == 0) cout << (index*5) << "%";
        else cout << ".";
        cout.flush();
        ++index;
    }
}

// Finishes the status bar
void finish_statusbar(int & index) {

    write_statusbar(index, 0, 1); // This will write all remaining percentages
    cout << "100%";
    cout.flush();
}

// Tries to find the specified value in the array and replaces it.
// Important: assumes that the value exists!
void find_and_replace(edge* array, edge find, edge replace){
    long i = 0;
    while (array[i]!=find)
		++i;
    array[i]= replace;
}


/*********************
 *  MAIN-Definitions *
 *********************/

// read a graph from "filename" and return it.
maingraph read_graph(char *filename, bool directed)
{
    maingraph ret;
    FILE *in = fopen(filename, "r");
    if (in != 0) {
	cout << " - Reading file \'" << filename << "\'." << endl;
	// Initialize variables
	const string numbers = "0123456789";
	ret.n = 0;
	ret.num_nodes = 0;
	ret.num_lonely_nodes = 0;
	bool error_flag = false;
	vertex u;
	vertex v;
	long weight;
	long linenumber = 0;
	edge e;
	int flag;
	char buffer[40];
	// Parse input File
	while (fgets(buffer, 100, in) != NULL) {
	    flag = sscanf(buffer, "%d %d %d", &u, &v, &weight);
	    if (flag > 1) {
			if (u > ret.n) {
				ret.n = u;
				if (ret.num_neighbours.size() <= ret.n) {
					ret.num_neighbours.resize(ret.n+1);
                    ret.num_undir.resize(ret.n+1);
				}
			}
			if (v > ret.n) {
				ret.n = v;
				if (ret.num_neighbours.size() <= ret.n) {
					ret.num_neighbours.resize(ret.n+1);
                    ret.num_undir.resize(ret.n+1);
				}
			}
			if (u != v) {
				if (u < v) {
					e = new_edge(u, v);
					if (ret.edges.find(e) == ret.edges.end()) {
						ret.num_neighbours[u]++;
						ret.num_neighbours[v]++;
						ret.edges[e] = DIR_U_T_V;
					} else {
						ret.edges[e] |= DIR_U_T_V;
					}
					if (!directed) {
						ret.edges[e] = UNDIR_U_V;
					}
				} else { // v < u
					e = new_edge(v, u);
					if (ret.edges.find(e) == ret.edges.end()) {
						ret.num_neighbours[u]++;
						ret.num_neighbours[v]++;
						ret.edges[e] = DIR_V_T_U;
					} else {
						ret.edges[e] |= DIR_V_T_U;
					}
					if (!directed) {
						ret.edges[e] = UNDIR_U_V;
					}
				}
                if (ret.edges[e] == UNDIR_U_V) {
                    ++ret.num_undir[u];
                    ++ret.num_undir[v];
                }
			}
	    } else {
			cerr << "ERROR: Wrong format for input file \'"
				<< filename << "\' in line " << linenumber
				<< "." << endl << "       ";
			throw domain_error("Expected line format is <int> <int> [<int>].");
		}
		++linenumber;
	}
    } else {
		string err = "ERROR: Unable to open file \'"
			         + string(filename) + "\' for input.";
		throw domain_error(err);
    }
    ++ret.n;
    fclose(in);
    ret.directed = directed;
    return ret;
}

// Allocates the arrays of the maingraph
void alloc_graph_arrays(maingraph & maing)
{
    maing.neighbours = new vertex *[maing.n+1];
	maing.num_larger_neighbours = new long[maing.n+1];
    maing.v_util = new unsigned long[maing.n+1];
    for (vertex i = 0; i != maing.n; ++i) {
      if (maing.directed)
          maing.neighbours[i] = new vertex[maing.num_neighbours[i]+maing.num_undir[i]];
      else
          maing.neighbours[i] = new vertex[maing.num_neighbours[i]];
    }
}

void build_graph(maingraph & maing)
{
    for (vertex i = 0; i != maing.n; ++i) {
      maing.v_util[i] = NILLVERTEX;          // reset util fields;
	  maing.num_larger_neighbours[i] = 0;
      maing.num_neighbours[i] = 0;           // num_neighbours might have changed during randomize
    }
    maing.m = maing.edges.size();       // Number of edges might have changed during "no regard"-randomization
	maing.maxnumlargerneighbors = 0;
    maing.num_undir_edges = 0;
    maing.num_dir_edges = 0;

    for (hash_map < edge, edgetype >::const_iterator iter = maing.edges.begin();
	 iter != maing.edges.end(); ++iter) {
		vertex u = edge_get_u(iter->first);
		vertex v = edge_get_v(iter->first);
		edgetype etype = iter->second;
        if (etype == UNDIR_U_V) {
			++maing.num_undir_edges;
        } else {
			++maing.num_dir_edges;
		}
		(maing.neighbours[u])[maing.num_neighbours[u]] = v;
		++maing.num_neighbours[u];
		(maing.neighbours[v])[maing.num_neighbours[v]] = u;
		++maing.num_neighbours[v];
		if (u>v) {
			++maing.num_larger_neighbours[v];
			if (maing.num_larger_neighbours[v]> maing.maxnumlargerneighbors)
			    ++maing.maxnumlargerneighbors;
		} else { // v>u
			++maing.num_larger_neighbours[u];
			if (maing.num_larger_neighbours[u]> maing.maxnumlargerneighbors)
			    ++maing.maxnumlargerneighbors;
		}
    }

	//Sort the neighbour array for every vertex
	for (int i = 0; i != maing.n; ++i)
		std::sort(maing.neighbours[i], maing.neighbours[i] + maing.num_neighbours[i]);

	return;
}

// estimates the size of the subgraph tree using TREESMPLS samples.
uint64 est_tree_size(const maingraph & g, long* v_extension, uint64 TREESMPLS, short G_N)
{
   	randlib::rand rand(time(NULL)); // Init Random-Generator
	register long* v_subgraph = new long[G_N+1];
	register long counter;
	register vertex SOURCE_VERTEX; 
	register vertex FILL_VERTEX;
	register long* min_scope = new long[G_N+1];
	register long* max_scope = new long[G_N+1];
	register long* scope_place = new long[G_N+1];
	register short depth;

    uint64 APPROX_TREE_SIZE = 0;
	for (uint64 th = 0; th!=TREESMPLS; ++th)  {
		vertex v = ++rand % g.n;
		uint64 approx = 1;
		bool flag = false;
		if (g.num_neighbours[v] > 0) {
			depth = 1;
			min_scope[depth] = 0;
			max_scope[depth] = 1;
			scope_place[depth] = 0;
			v_extension[0] = v;
			v_subgraph[depth-1] = NILLVERTEX;
			g.v_util[v] = v;
			while (depth != 0)  {
				if (flag) { //scope_place[depth] == max_scope[depth] || flag)	{
					if (min_scope[depth] != max_scope[depth]) {
						counter = max_scope[depth] - 1;
						while (counter >= 0 && g.v_util[v_extension[counter]] == v_subgraph[depth-1]) {				
							g.v_util[v_extension[counter]] = NILLVERTEX;
							--counter;
						}
					}
					--depth;
				} else {
					if (depth == G_N) { 
						flag = true; 
						while (scope_place[depth] != max_scope[depth]) {
							v_subgraph[depth] = v_extension[scope_place[depth]];
							++scope_place[depth];
						}
					} else {  
						if (scope_place[depth] == max_scope[depth]) {
							approx = 0;
							flag = true;
						} else {
							scope_place[depth] += (++rand) % (max_scope[depth]-scope_place[depth]); 
							SOURCE_VERTEX = v_extension[scope_place[depth]];
							++scope_place[depth];
							v_subgraph[depth] = SOURCE_VERTEX;
							min_scope[depth+1] = max_scope[depth];
							scope_place[depth+1] = scope_place[depth];
							max_scope[depth+1] = min_scope[depth+1];
							++depth;
							counter = g.num_neighbours[SOURCE_VERTEX] - 1;
							if (counter > 0)   {
								FILL_VERTEX = g.neighbours[SOURCE_VERTEX][counter];
							}
							while (counter > -1 && g.neighbours[SOURCE_VERTEX][counter] > v)	{
								FILL_VERTEX = g.neighbours[SOURCE_VERTEX][counter];
								if (g.v_util[FILL_VERTEX] == NILLVERTEX)	{
									v_extension[max_scope[depth]] = FILL_VERTEX;
									g.v_util[FILL_VERTEX] = SOURCE_VERTEX;
									++max_scope[depth];
								}
								--counter;
							}
							approx *= (max_scope[depth] - scope_place[depth]);
						}
					}
				}
			}
			g.v_util[v] = NILLVERTEX;
			APPROX_TREE_SIZE += approx;
		}

	}
    delete[] v_subgraph;
    delete[] min_scope;
    delete[] max_scope;
    delete[] scope_place;

    return APPROX_TREE_SIZE;
}

// Samples / Enumerates the graph returns results in result_graphs, and the time directly.
double sampling(const maingraph & maing, long* v_extension, short G_N, const int G_M,
                bool fullenumeration, const double* prob, const uint64 equiv100p,
                const bool write_status, graph* canon, graph* g, set* gv, int* lab,
                int* ptn, int* orbits, optionblk & options, statsblk & stats,
                setword* workspace, hash_map < graph64, uint64 > & result_graphs,
                uint64 & count_subgr)
{
    // Init for Statusbar
    if (write_status){
        cout << "0%";
        cout.flush();
    }
    register int percentage_index = 1;

    // Init for main loop
	count_subgr = 0;
	randlib::rand rand(time(NULL));
	register long* v_subgraph = new long[G_N+1];
	register long counter;
	register vertex SOURCE_VERTEX; 
	register vertex FILL_VERTEX;
	register long* min_scope = new long[G_N+1];
	register long* max_scope = new long[G_N+1];
	register long* scope_place = new long[G_N+1];
	register short depth;

  	register vertex* start_vertices = new vertex[maing.n+2];//##%!
	register short* scope_place_point = new short[G_N+1];//##%!
	register unsigned long** scope_place_loc = new unsigned long*[G_N+1];

    // Beginning actual sampling
    for (int i = 0; i != G_N+1; ++i) {
		scope_place_loc[i] = new unsigned long[maing.maxnumlargerneighbors*G_N];
	}
	//scope_place_loc[8][maxnumlargerneighbors*G_N];//##%!
	//gen_selection(0,tst,0.174,arr,rand);
	if (fullenumeration) {
		for (vertex v = 0; v != maing.n; ++v)
			start_vertices[v] = v;
		start_vertices[maing.n] = maing.n;
	} else {
		gen_selection(0,maing.n,prob[0],start_vertices,rand);
	}
	
	
	clock_t start_time(clock());
	for (int th = 0; th!=1; ++th) // this line for time measurements only 
	{

	int idx = 0;
	while (start_vertices[idx]!=maing.n)
	{	
		vertex v = start_vertices[idx];
		++idx;
		//if (fullenumeration || ((++rand)%MAXPROB) <= prob[0])///////////////////////////////////////////////
		{
		
		depth = 1;
		min_scope[depth] = 0;
		max_scope[depth] = 1;
		scope_place[depth] = 0;
		v_extension[0] = v;
		v_subgraph[depth-1] = NILLVERTEX;
		maing.v_util[v] = v;
		scope_place_point[depth] = 0;
		scope_place_loc[1][0] = 0;
		scope_place_loc[1][1] = 1;

		while (depth != 0)  {
			if (scope_place[depth] == max_scope[depth])	{ // go lower
				if (min_scope[depth] != max_scope[depth]) {
					counter = max_scope[depth] - 1;
					while (counter >= 0 && maing.v_util[v_extension[counter]] == v_subgraph[depth-1]) {
						maing.v_util[v_extension[counter]] = NILLVERTEX;
						--counter;
					}
				}
				--depth;
			} else {
				if (depth == G_N) { 
					while (scope_place[depth] != max_scope[depth])
					{
						{
						++count_subgr;
						v_subgraph[depth] = v_extension[scope_place[depth]];
						// edit nauty-graph 
						for (int i = 0; i != depth-1; ++i) {
							vertex uc = v_subgraph[i+1];
							vertex vc = v_subgraph[depth];
							edge e_check = edge_code(uc, vc);
							DELELEMENT(GRAPHROW(g, i, G_M), depth-1);
							DELELEMENT(GRAPHROW(g, depth-1, G_M), i);
							if (maing.edges.find(e_check) != maing.edges.end()) {
								edgetype ec = maing.edges.find(e_check)->second;
								if (uc < vc) {
									switch (ec) {
										case DIR_U_T_V:
										gv = GRAPHROW(g, i, G_M);
										ADDELEMENT(gv, depth-1);
										break;
										case DIR_V_T_U:
										gv = GRAPHROW(g, depth-1, G_M);
										ADDELEMENT(gv, i);
										break;
										case UNDIR_U_V:
										gv = GRAPHROW(g, i, G_M);
										ADDELEMENT(gv, depth-1);
										gv = GRAPHROW(g, depth-1, G_M);
										ADDELEMENT(gv, i);
										break;
									}
								} else {
									switch (ec) {
										case DIR_U_T_V:
										gv = GRAPHROW(g, depth-1, G_M);
										ADDELEMENT(gv, i);
										break;
										case DIR_V_T_U:
										gv = GRAPHROW(g, i, G_M);
										ADDELEMENT(gv, depth-1);
										break;
										case UNDIR_U_V:
										gv = GRAPHROW(g, depth-1, G_M);
										ADDELEMENT(gv, i);
										gv = GRAPHROW(g, i, G_M);
										ADDELEMENT(gv, depth-1);
										break;
									}
								}
							}
							
						}
						// end edit nauty-graph

						// Canonical labelling with nauty
						nauty(g, lab, ptn, NILSET, orbits, &options,
						  &stats, workspace, 160 * MAXM, G_M, G_N,
						  canon);
						graph64 res_gr = 0ULL;
						for (int a = 0; a < G_N; ++a) {
							gv = GRAPHROW(canon, a, G_M);
							for (int b = 0; b < G_N; ++b) {
								res_gr <<= 1;
								if (ISELEMENT(gv, b)) {
									res_gr |= 1;
								}
							}
						}
						result_graphs[res_gr]++;
						// end nauty
						if (fullenumeration)
						{
							++scope_place[depth];
						} else {
							++scope_place_point[depth];
							scope_place[depth] = scope_place_loc[depth][scope_place_point[depth]];						
						}
						}

					}

					// Progress indicator
                    if (write_status) write_statusbar(percentage_index, equiv100p, count_subgr);

				} else {  //go deeper
					SOURCE_VERTEX = v_extension[scope_place[depth]];
					if (fullenumeration) {
						++scope_place[depth];
					} else {
						++scope_place_point[depth];
						scope_place[depth] = scope_place_loc[depth][scope_place_point[depth]];
					}

					{//!!

					v_subgraph[depth] = SOURCE_VERTEX;
					// edit nauty-graph 
					for (int i = 0; i != depth-1; ++i) {
						vertex uc = v_subgraph[i+1];
						vertex vc = SOURCE_VERTEX;
						edge e_check = edge_code(uc, vc);
						DELELEMENT(GRAPHROW(g, i, G_M), depth-1);
						DELELEMENT(GRAPHROW(g, depth-1, G_M), i);
						if (maing.edges.find(e_check) != maing.edges.end()) {
							edgetype ec = maing.edges.find(e_check)->second;
							if (uc < vc) {
								switch (ec) {
									case DIR_U_T_V:
									gv = GRAPHROW(g, i, G_M);
									ADDELEMENT(gv, depth-1);
									break;
									case DIR_V_T_U:
									gv = GRAPHROW(g, depth-1, G_M);
									ADDELEMENT(gv, i);
									break;
									case UNDIR_U_V:
									gv = GRAPHROW(g, i, G_M);
									ADDELEMENT(gv, depth-1);
									gv = GRAPHROW(g, depth-1, G_M);
									ADDELEMENT(gv, i);
									break;
								}
							} else {
								switch (ec) {
									case DIR_U_T_V:
									gv = GRAPHROW(g, depth-1, G_M);
									ADDELEMENT(gv, i);
									break;
									case DIR_V_T_U:
									gv = GRAPHROW(g, i, G_M);
									ADDELEMENT(gv, depth-1);
									break;
									case UNDIR_U_V:
									gv = GRAPHROW(g, depth-1, G_M);
									ADDELEMENT(gv, i);
									gv = GRAPHROW(g, i, G_M);
									ADDELEMENT(gv, depth-1);
									break;
								}
							}
						}
					}
					// end edit nauty-graph 
					min_scope[depth+1] = max_scope[depth];
					scope_place[depth+1] = scope_place[depth];
					max_scope[depth+1] = min_scope[depth+1];
					++depth;
					counter = maing.num_neighbours[SOURCE_VERTEX] - 1;
					if (counter > 0)   {
						FILL_VERTEX = maing.neighbours[SOURCE_VERTEX][counter];
					}
					while (counter > -1 && maing.neighbours[SOURCE_VERTEX][counter] > v)	{
						FILL_VERTEX = maing.neighbours[SOURCE_VERTEX][counter];
						if (maing.v_util[FILL_VERTEX] == NILLVERTEX)	{
							v_extension[max_scope[depth]] = FILL_VERTEX;
							maing.v_util[FILL_VERTEX] = SOURCE_VERTEX;
							++max_scope[depth];
						}
						--counter;
					}
					}//!!
					
					//GENERATE SELECTION AND RESET SCOPE PLACE POINT
					if (! fullenumeration)
					{
						scope_place_point[depth] = 0;
						gen_selection(scope_place[depth], max_scope[depth],
						              prob[depth-1],scope_place_loc[depth],rand);
						scope_place[depth] = scope_place_loc[depth][0];
					}
				}
			}
		}
		maing.v_util[v] = NILLVERTEX;
		}
	}

	}

    // Finish the statusbar
    if (write_status) finish_statusbar(percentage_index);

    delete[] v_subgraph;
    delete[] min_scope;
    delete[] max_scope;
    delete[] scope_place;
    delete[] start_vertices;
    delete[] scope_place_point;
	for (int i = 0; i != G_N+1; ++i) {
		delete[] scope_place_loc[i];
	}
    delete[] scope_place_loc;

   	// Done with sampling / enumeration: stop the clock, return the time
    return double (clock() - start_time) / CLOCKS_PER_SEC;
}



double randomize_graph(maingraph & maing, short random_type,
                     int num_exchanges, int num_tries, long & total_tries,
                     long & total_success)
{
    randlib::rand rand(time(NULL)); // Init Random-Generator
    edge* edgelist; // The edges are written into this array
    long* open_positions; // Stores the free places in the edgelist array.
    long op_index = 0;    // The index one behind the last open position.

    edge e, p, ne1, ne2, ne3, single, undir, third;
    // e is the edge, which needs an exchange partner, p is this (random) found partner,
    // ne1-3 are newly created edges (often it is first checked, whether those edges already exist)
    // single, undir and third are used for the global complex rules.
    // you start with an single/directed edge, pick an undirected edge, and a third directed edge.
    vertex eu, ev, pu, pv, unu, unv, tu, tv, begin_path, middle, end_path, other;
    // the first are the u and v of the above defined vertices.
    // begin_path, middle, end_path and other are again used for the complex rules.
    // there is an path between three vertices, and the fourth is the "other" vertex.
    edgetype et, pt, est, pst, ne1t, ne2t, st, tt;
    // Types of the edges above: st is type of edge "single", tt type of the edge "third".
    // est and pst are used for "no regard"-graphs: edge split type and partner split type
    hash_map < edge, edgetype >::const_iterator find_ne1, find_ne2;
    // those store the values of finding ne1 and ne2 in the hashmap.
    bool success, found_path;
    // success is set if an exchange was made, found_path, if an path between the vertices is found
    // These longs store the numbers the edges have in the edgelist array.
    long p_number, undir_number, single_number;
    // These vars do some statistics
    long complex4used = 0, complex4tried = 0,
         complex5used = 0, complex5tried = 0, failedbidir = 0;

    // In NO_REGARD-Graphs the number of edges can increase due to an edge-split
    if (random_type == NO_REGARD){
        edgelist = new edge[maing.m+maing.num_undir_edges];
        open_positions = new long[maing.m];
    } else
        edgelist = new edge[maing.m];

	clock_t start_time(clock()); // Start the time measurement

    // Turn edges hashmap into an array. It is sorted: First dir edges, then undir edges
    { long i=0, j=maing.num_dir_edges;
    for (hash_map < edge, edgetype >::const_iterator iter = maing.edges.begin();
         iter != maing.edges.end(); ++iter) {
        if (iter->second == UNDIR_U_V){
            edgelist[j] = iter->first;
            ++j;
        } else {
            edgelist[i] = iter->first;
            ++i;
        }
    }
    }

    for (int exchange_ctr = 0; exchange_ctr < num_exchanges; ++exchange_ctr)
    {
        // for each edge: try to find a partner
        for (long e_number = 0; e_number < maing.m; ++e_number){
            e = edgelist[e_number];  // select edge according to counter
            if (e == 0xFFFFFFFFFFFFFFFFULL) { // if the edge has been erased (which can happen in NO_REGARD-Graphs only)
                continue;                     // continue the for-loop at the next edge.
            }
            eu = edge_get_u(e);  // do the lookups for e
            ev = edge_get_v(e);
            et = maing.edges[e];
            success = false;       // reset success flag
            // Try to use a rule until one is successful or num_tries is exceeded
            for (int tries_ctr = 0; (!success) && (tries_ctr < num_tries); ++tries_ctr){
                found_path = false;  // reset found_path flag
                ++total_tries;       // count tries for the statistic
                switch (random_type){
// NO_REGARD ------------
                case NO_REGARD:
                    p_number = ++rand % maing.m;
                    p = edgelist[p_number]; // find partner edge
                    while (p == 0xFFFFFFFFFFFFFFFFULL) { // if the edge has been erased
                        p_number = ++rand % maing.m;     // chose a new one.
                        p = edgelist[p_number];
                    }
                    pu = edge_get_u(p);
                    pv = edge_get_v(p);
                    pt = maing.edges[p];
                    if (pu != eu && pv != ev && pu != ev && pv != eu){ // ensure that vertices are disjunct
                        // In case e or p are bidir. choose only one direction.
                        est = (et == UNDIR_U_V) ? (++rand % 2)+1 : et; // split undirected edge
                        pst = (pt == UNDIR_U_V) ? (++rand % 2)+1 : pt;
                        if (est == pst) { // Determine directions and type of new edges
                            if (eu < pv){ // This creates ne1 as an new edge from eu to pv
                                ne1 = new_edge(eu, pv);
                                ne1t = est;
                            } else {      // In case pv is smaller than eu the direction has to be reversed.
                                ne1 = new_edge(pv, eu);
                                ne1t = reverse(est);
                            }
                            if (pu < ev){ // ne2 and its type are defined going from pu to ev.
                                ne2 = new_edge(pu, ev);
                                ne2t = pst;
                            } else {
                                ne2 = new_edge(ev, pu);
                                ne2t = reverse(pst);
                            }
                        } else { // est != pst In case the types of e and p are unequal,
                                 // p is mentally turned, so that the new edges now
                                 // lead from eu to pu and from pv to ev.
                            if (eu < pu){
                                ne1 = new_edge(eu, pu);
                                ne1t = est;
                            } else {
                                ne1 = new_edge(pu, eu);
                                ne1t = reverse(est);
                            }
                            if (pv < ev){
                                ne2 = new_edge(pv, ev);
                                ne2t = reverse(pst); // because p is handled as if reversed, we have to reverse the type of the new edge.
                            } else {
                                ne2 = new_edge(ev, pv);
                                ne2t = pst; // This is really reverse(reverse(pst))!
                            }
                        } // end if types equal
                        find_ne1 = maing.edges.find(ne1); // Look up the edges we want to create.
                        find_ne2 = maing.edges.find(ne2);
                        // Check if they are not already existing.
                        if ((find_ne1 == maing.edges.end() || find_ne1->second == reverse(ne1t)) &&
                            (find_ne2 == maing.edges.end() || find_ne2->second == reverse(ne2t))){
                            // Add to hashmap and update edgelist.
                            // The types are not regarded, as the list does not need to be sorted for this random-type
                            if (maing.edges.find(ne1) == maing.edges.end()){ // ne1 does not exist:
                                maing.edges[ne1] = ne1t;         // Create entry for ne1
                                if (et == est){                  // In this case, e has to be erased:
                                    edgelist[e_number] = ne1;    // Overwrite e with ne1
                                } else {                         // in this case, e should stay in the list:
                                    if (op_index == 0){          // There are no open positions:
                                        edgelist[maing.m] = ne1; // Add ne1 at the end of list.
                                        ++maing.m;
                                    } else {                     // There are open positions:
                                        --op_index;              // Write ne2 in one of them.
                                        edgelist[open_positions[op_index]] = ne1;
                                    }
                                }
                            }
                            else {                               // ne1 exists:
                                maing.edges[ne1] = UNDIR_U_V;    // ne1 is now undirected.
                                if (et == est){                  // In this case, e has to be erased:
                                    edgelist[e_number] = 0xFFFFFFFFFFFFFFFFULL; // This is an kind of erase...
                                    open_positions[op_index] = e_number; // The position of e is now open.
                                    ++op_index;
                                }
                                // If (et != est) e should remain in the list and ne1 is already in the list:
                                // There is no work at all.
                            }

                            if (maing.edges.find(ne2) == maing.edges.end()){ // ne2 does not exist:
                                maing.edges[ne2] = ne2t;         // Create entry for ne2
                                if (pt == pst){                  // In this case, p has to be erased:
                                    edgelist[p_number] = ne2;    // Overwrite p with ne2
                                } else {                         // in this case, p should stay in the list:
                                    if (op_index == 0){          // There are no open positions:
                                        edgelist[maing.m] = ne2; // Add ne2 at the end of list.
                                        ++maing.m;
                                    } else {                     // There are open positions:
                                        --op_index;              // Write ne2 in one of them.
                                        edgelist[open_positions[op_index]] = ne2;
                                    }
                                }
                            }
                            else {                               // ne2 exists:
                                maing.edges[ne2] = UNDIR_U_V;    // ne2 is now undirected.
                                if (pt == pst){                  // In this case, p has to be erased:
                                    edgelist[p_number] = 0xFFFFFFFFFFFFFFFFULL; // This is an kind of erase...
                                    open_positions[op_index] = p_number; // The position of p is now open.
                                    ++op_index;
                                }
                                // If (pt != pst) p should remain in the list and ne2 is already in the list:
                                // There is no work at all.
                            }

                            // Erase from Hashmap:
                            if (et == est)                    // In this case e has to be erased.
                                maing.edges.erase(e);
                            else {                            // et != est: e has to stay in the map
                                maing.edges[e] = reverse(est);// Set e only to the not switched dir.
                            }

                            if (pt == pst)                    // In this case p has to be erased.
                                maing.edges.erase(p);         // pt != pst: p has to stay in the map
                            else {
                                maing.edges[p] = reverse(pst);// Set p only to the not switched dir.
                            }
                            success=true;
                        }  // end if edges do not exist
                    } // end if vertices unequal
                break;
                case GLOBAL_CONST:
// GLOBAL_CONST
// Complex 4 vertex rule
                    // If there are undir edges and directed edges, with an 50% probability:
                    // Use complex 4-vertex-rule (probability can be adjusted)
                    if (maing.num_undir_edges != 0 && maing.num_dir_edges != 0
                        && ++rand % 2 == 0){
                        ++complex4tried; // Count tries of this rule
                        // If edge e is undirected: randomly select directed edge
                        if (et == UNDIR_U_V){
                            single_number = ++rand % maing.num_dir_edges;
                            undir_number = e_number; // store the numbers to be able to update edgelist array
                            single = edgelist[single_number];
                            st = maing.edges[single]; // lookup the single edge
                            // define the middle and the end of the path according to the direction of "single"
                            middle = (st == DIR_U_T_V) ? edge_get_u(single) : edge_get_v(single);
                            end_path = (st == DIR_U_T_V) ? edge_get_v(single) : edge_get_u(single);
                            undir = e; // edge e is from now on called "undir"
                            unu = eu;  // its vertices are unu and unv
                            unv = ev;
                        } else {  // If edge e is directed randomly select undirected edge
                            undir_number = (++rand % maing.num_undir_edges) + maing.num_dir_edges;
                            single_number = e_number; // store the numbers to be able to update edgelist array
                            undir = edgelist[undir_number]; // the undirected edge is looked up
                            unu = edge_get_u(undir);        // and its vertices defined
                            unv = edge_get_v(undir);
                            single = e;                     // edge e is from now on called "single"
                            st = et;                        // its type is st (It has to be one of the dir. types)
                            middle = (st == DIR_U_T_V) ? eu : ev;    // middle and end_path are defined
                            end_path =  (st == DIR_U_T_V) ? ev : eu; // according to the direction of edge "single"
                        }
                        // Ensure that no two vertices are the same:
                        if (unu != middle && unv != end_path && unu != end_path && unv != middle){
                            begin_path = unu;   // first try: unu is the begin of the path
                            other = unv;        // unv is the "other" vertex, which is not in the path
                            // Define u and v of the "third" edge
                            if (unu < middle) {
                               tu = begin_path;
                               tv = middle;
                            } else {
                               tu = middle;
                               tv = begin_path;
                            }
                            third = new_edge(tu, tv); // Define the third edge.

                            if (maing.edges.find(third) != maing.edges.end()) { // third should exist
                                tt = maing.edges[third];
                                if (tt == DIR_U_T_V && tu == begin_path) found_path = true; // third is pointing into the right direction
                                if (tt == DIR_V_T_U && tv == begin_path) found_path = true; // a path is found.
                            }
                            if (!found_path){ // If no path could be found with unu as begin path:
                                begin_path = unv; // Use unv as begin
                                other = unu;      // unu then is the other vertex, not in the path
                                if (unv < middle) {  // same as above
                                    tu = begin_path;
                                    tv = middle;
                                } else {
                                    tu = middle;
                                    tv = begin_path;
                                }
                                third = new_edge(tu, tv); // third is looked up again
                                if (maing.edges.find(third) != maing.edges.end()) {
                                   tt = maing.edges[third];  // as above: Is third pointing into the right direction?
                                   if (tt == DIR_U_T_V && tu == begin_path) found_path = true;
                                   if (tt == DIR_V_T_U && tv == begin_path) found_path = true;
                               }
                            } // end if use unv as begin
                            if (found_path) { // if an path was found:
                                if (other < end_path){ // Define ne1 as pointing from other to end path.
                                    ne1 = new_edge(other, end_path);
                                    ne1t = DIR_U_T_V;
                                 } else {
                                    ne1 = new_edge(end_path, other);
                                    ne1t = DIR_V_T_U;
                                }
                                // if this new edge ne1 does not exist: begin the switching:
                                if (maing.edges.find(ne1) == maing.edges.end()){
                                    maing.edges[ne1] = ne1t;        // ne1 is created
                                    maing.edges[third] = UNDIR_U_V; // the third edge is nor undirected
                                    if (begin_path < other){  // The undirected edge is now pointing from begin_path to other
                                        maing.edges[undir] = DIR_U_T_V;
                                    } else {
                                        maing.edges[undir] = DIR_V_T_U;
                                    }

                                    maing.edges.erase(single); // the single edge is erased.
                                    // ne1 is an dir edge, therefore it can take the place of single in the edgelist
                                    edgelist[single_number]= ne1;
                                    // undir has become directed, therefore it
                                    // takes the place of third, which used to be directed
                                    find_and_replace(edgelist, third, undir);
                                    // third is now undirected and placed at "undir"s position.
                                    edgelist[undir_number] = third;
                                    success = true;
                                    ++complex4used; // The rule was used successfully
                                } // end if switching
                            } // end if found path
                        } // end if vertices unequal

// Easy rules -------------
                    } else { // Use easy rules or 5-vertex-rule

                    // If edge e is undirected: randomly select undirected edge
                    if (et == UNDIR_U_V){ // randomly select undirected edge
                       p_number = (++rand % maing.num_undir_edges) + maing.num_dir_edges;
                    } else { // et = DIR: select another directed edge
                       p_number = ++rand % maing.num_dir_edges;
                    } // end if et = UNDIR
                    p = edgelist[p_number]; // find partner edge
                    pu = edge_get_u(p);
                    pv = edge_get_v(p);
                    pt = maing.edges[p];    // lookup type of partner edge
                    if (pu != eu && pv != ev && pu != ev && pv != eu){
                        if (et == pt) { // Determine directions and type of new edges
                            if (eu < pv){  // These are the same procedures as in "NO REGARD" above
                                ne1 = new_edge(eu, pv);  // New edge ne1 from eu to pv
                                ne1t = et;
                            } else {
                                ne1 = new_edge(pv, eu);
                                ne1t = reverse(et);
                            }
                            if (pu < ev){                // New edge ne2 from pu to ev
                                ne2 = new_edge(pu, ev);
                                ne2t = pt;
                            } else {
                                ne2 = new_edge(ev, pu);
                                ne2t = reverse(pt);
                            }
                        } else { // et != pt Two directed edges point in opposite directions.
                        // as above p is then mentally turned.
                            if (et==UNDIR_U_V)  // Assertion: This should never ever happen.
                                cout << "Tried to change dir with undir edge!! Big trouble" << endl;
                            if (eu < pu){               // New edge ne1 from eu to pu
                                ne1 = new_edge(eu, pu);
                                ne1t = et;
                            } else {
                                ne1 = new_edge(pu, eu);
                                ne1t = reverse(et);
                            }
                            if (pv < ev){               // New edge ne2 from pv to ev
                                ne2 = new_edge(pv, ev);
                                ne2t = reverse(pt);
                            } else {
                                ne2 = new_edge(ev, pv);
                                ne2t = pt;
                            }
                        } // end if types equal

                        // Look up the edges we want to create
                        find_ne1 = maing.edges.find(ne1);
                        find_ne2 = maing.edges.find(ne2);
                        // Check whether they are not already existing:
                        if (find_ne1 == maing.edges.end() &&
                            find_ne2 == maing.edges.end()) {
                            maing.edges.erase(e);   //Update Hashmap
                            maing.edges.erase(p);   // e and p are erased.
                            maing.edges[ne1] = ne1t; // ne1 and ne2 added with their types
                            maing.edges[ne2] = ne2t;
                            edgelist[e_number]= ne1;         // Update edgelist
                            edgelist[p_number] = ne2; // ne1 takes the position of e, ne2 the position of p
                            success=true;
// 5 Vertex Rule ---------------
                        // There is an connection between the vertices: try to use 5-vertex-rule
                        // Edge e has to be dir as it shall be part of the path.
                        // There has to be at least one undir edge and at least two dir edges, to use the rule
                        } else if (et != UNDIR_U_V && maing.num_undir_edges != 0 && maing.num_dir_edges > 1) {
                            ++complex5tried; // Count tries of this rule
                            if (find_ne1 != maing.edges.end()) { // If edge ne1 exists:
                                // try to use it as the third edge for the switching rule
                                // The name third comes from the 4 vertex rule, where this is really the third edge to be picked.
                                third = ne1;
                                // edges changed are: e, undir (defined later), and third
                                tt = find_ne1->second; // the type of third according to the hashtable.
                                if (et == DIR_U_T_V){
                                   // if edge e points from u to v, it is the second edge in the path,
                                   // because ne1 has eu as one vertex, so eu is the middle vertex of the path.
                                   if (edge_get_u(third) == eu){ // If third and e have the same "u"
                                       if (tt == DIR_V_T_U){ // and the directions are ok
                                           found_path = true;                // we have found an path
                                           begin_path = edge_get_v(third);   // begin, middle and end are set accordingly
                                           middle = eu;
                                           end_path = ev;

                                       }
                                   } else { // edge_get_v(third) == eu
                                       // third and e do not have the same u =>
                                       // third has to point in the opposite direction for the path to be correct.
                                       if (tt == DIR_U_T_V){
                                           found_path = true;
                                           begin_path = edge_get_u(third);
                                           middle = eu;
                                           end_path = ev;
                                       }
                                   }
                                } else { // et == DIR_V_T_U
                                   // if edge e points from v to u, it is the first edge in the path.
                                   if (edge_get_u(ne1) == eu){
                                       if (find_ne1->second == DIR_U_T_V){
                                           found_path = true;
                                           begin_path = ev;
                                           middle = eu;
                                           end_path = edge_get_v(third);
                                       }
                                   } else { // edge_get_v(third) == eu
                                       if (find_ne1->second == DIR_V_T_U){
                                           found_path = true;
                                           begin_path = ev;
                                           middle = eu;
                                           end_path = edge_get_u(third);
                                       }
                                   }
                                } // end if lookup of et
                            } // end if "ne1 exists"
                            if (!found_path && find_ne2 != maing.edges.end()) {
                            // If we haven't already found a path and ne2 exists,
                            // we are trying to use it as "third" edge.
                                third = ne2; // edges changed are: e, undir (defined later), and third
                                // ne2 certainly has ev as one vertex.
                                tt = find_ne2->second; // Type of third according to the hashtable
                                if (et == DIR_U_T_V){  // These things are going the same way, as with ne1 above
                                   if (edge_get_u(third) == ev){
                                       if (tt == DIR_U_T_V){
                                           found_path = true;
                                           begin_path = eu;
                                           middle = ev;
                                           end_path = edge_get_v(third);
                                       }
                                   } else { // edge_get_v(third) == eu
                                       if (tt == DIR_V_T_U){
                                           found_path = true;
                                           begin_path = eu;
                                           middle = ev;
                                           end_path = edge_get_u(third);
                                       }
                                   }
                                } else { // et == DIR_V_T_U
                                   if (edge_get_u(third) == ev){
                                       if (tt == DIR_V_T_U){
                                           found_path = true;
                                           begin_path = edge_get_v(third);
                                           middle = ev;
                                           end_path = eu;
                                       }
                                   } else { // edge_get_v(ne2) == ev
                                       if (tt == DIR_U_T_V){
                                           found_path = true;
                                           begin_path = edge_get_u(third);
                                           middle = ev;
                                           end_path = eu;
                                       }
                                   }
                                } // end if lookup of et
                            } // end if ne2 exists
                        } // end if new edges do not exist
                     } // end if vertices unequal (These brackets finish the easy rule usage)
                     if (found_path){ // Path was found in section above: Using 5-vertex-rule
                     // Try to find an bidir edge fitting to the path num_tries times.
                     for (int tries_ctr2 = 0; (!success) && (tries_ctr2 < num_tries); ++tries_ctr2){
                         // randomly select an undir edge, number is stored once again for the edgelist update.
                         undir_number = (++rand % maing.num_undir_edges) + maing.num_dir_edges;
                         undir = edgelist[undir_number];
                         unu = edge_get_u(undir);
                         unv = edge_get_v(undir);
                         // if vertices are unequal
                         if (unu != begin_path && unu != middle && unu != end_path &&
                             unv != begin_path && unv != middle && unv != end_path) {
                             if (begin_path < unv) { // The new edge ne1 is pointing from begin_path to unv.
                                 ne1 = new_edge(begin_path, unv);
                                 ne1t = DIR_U_T_V;
                             } else {
                                 ne1 = new_edge(unv, begin_path);
                                 ne1t = DIR_V_T_U;
                             }
                             if (unv < end_path) {  // The new edge ne2 is pointing from unv to end_path.
                                 ne2 = new_edge(unv, end_path);
                                 ne2t = DIR_U_T_V;
                             } else {
                                 ne2 = new_edge(end_path, unv);
                                 ne2t = DIR_V_T_U;
                             }
                             if (unu < middle) {   // The new edge ne3 is undirected between unu and middle.
                                 ne3 = new_edge(unu, middle);
                             } else {
                                 ne3 = new_edge(middle, unu);
                             }
                             // Those three edges we want to create should not exist:
                             if (maing.edges.find(ne1) == maing.edges.end() &&
                                 maing.edges.find(ne2) == maing.edges.end() &&
                                 maing.edges.find(ne3) == maing.edges.end()) {
                                 // start the 5-vertice switching
                                 maing.edges[ne1] = ne1t;      // Insert ne1-3 with their types
                                 maing.edges[ne2] = ne2t;
                                 maing.edges[ne3] = UNDIR_U_V;
                                 maing.edges.erase(e);         // Erase the edges we worked with:
                                 maing.edges.erase(undir);     // e, undir and third
                                 maing.edges.erase(third);
                                 // third was undirected, ne1 is undir, so it can replace third.
                                 find_and_replace(edgelist, third, ne1);
                                 // e was undirected, ne2 ist undir, so it can replace e.
                                 edgelist[e_number] = ne2;
                                 // ne3 is undirected, so it can replace undir.
                                 edgelist[undir_number] = ne3;
                                 success = true;
                             } // end if new edges do not exist
                         } // end if vertices unequal
                         } // end for tries to find bidir edge
                         if (success) ++complex5used; // count the usage of this rule.
                         // if we were here, but did not succeed, we had an path, but no bidir edge for it.
                         else ++failedbidir;
                     } // end if found_path
                     } // end if which rule to use
                break;
// LOCAL_CONST
                case LOCAL_CONST: // The local const rules are exactly the same as the easy global const rules.
                                  // so comments become rare here.
                    if (et == UNDIR_U_V){ // randomly select undirected edge
                       p_number = (++rand % maing.num_undir_edges) + maing.num_dir_edges;
                    } else { // et = DIR
                       p_number = ++rand % maing.num_dir_edges;
                    } // end if et = UNDIR
                    p = edgelist[p_number]; // find partner edge
                    pu = edge_get_u(p);
                    pv = edge_get_v(p);
                    pt = maing.edges[p];
                    if (pu != eu && pv != ev && pu != ev && pv != eu){
                        if (et == pt) { // Determine directions and type of new edges
                            if (eu < pv){
                                ne1 = new_edge(eu, pv);
                                ne1t = et;
                            } else {
                                ne1 = new_edge(pv, eu);
                                ne1t = reverse(et);
                            }
                            if (pu < ev){
                                ne2 = new_edge(pu, ev);
                                ne2t = pt;
                            } else {
                                ne2 = new_edge(ev, pu);
                                ne2t = reverse(pt);
                            }
                        } else { // et != pt
                            if (et==UNDIR_U_V)  // This should never ever happen.
                                cout << "Tried to exchange dir with undir edge!!" << endl;
                            if (eu < pu){
                                ne1 = new_edge(eu, pu);
                                ne1t = et;
                            } else {
                                ne1 = new_edge(pu, eu);
                                ne1t = reverse(et);
                            }
                            if (pv < ev){
                                ne2 = new_edge(pv, ev);
                                ne2t = reverse(pt);
                            } else {
                                ne2 = new_edge(ev, pv);
                                ne2t = pt;
                            }
                        } // end if types equal
                        if (maing.edges.find(ne1)== maing.edges.end() &&
                            maing.edges.find(ne2)== maing.edges.end()) {
                            maing.edges.erase(e);   //Update Hashmap
                            maing.edges.erase(p);
                            maing.edges[ne1] = ne1t;
                            maing.edges[ne2] = ne2t;
                            edgelist[e_number]= ne1;         // Update edgelist
                            edgelist[p_number] = ne2;
                            success=true;
                        } // end if new edges do not exist
                     } // end if vertices unequal
                break;
                } // end switch
            } // end for tries
            if (success) total_success++; // count the amount of successful exchanges
        } // end for edges
    } // end for exchanges

    delete[] edgelist;
    if (random_type == NO_REGARD)
     delete[] open_positions;

    //Output Statistics
    //(remark: as total_tries and total_success are not initialised with 0, these are no longer correct)
    /* cout << "Randomize graph: " << endl;
    cout << "used " << total_tries << " tries to do " << total_success
         << " changes." << endl;
    if (random_type == GLOBAL_CONST) {
    cout << "tried to use complex 4 rule " << complex4tried << " times, "
         << complex4used << " were successful." << endl;
    cout << "tried to use complex 5 rule " << complex5tried << " times, "
         << complex5used << " were successful, " << endl;
    cout << failedbidir
         << " failed, because they did not find a proper bidir edge." << endl;
    }
    cout << endl; */

    build_graph(maing); // Update neighbour arrays
    return double (clock() - start_time) / CLOCKS_PER_SEC;

}
