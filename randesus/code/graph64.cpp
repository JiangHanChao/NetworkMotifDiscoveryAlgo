#include "graph64.hpp"

/*
 * Output Adjacency Matrix of 64 Bit Graph
 */
void adj(graph64 g) {
	short shift;
	for (int i= 0; i!=8 ; ++i) {
		for (int j= 0; j!=8 ; ++j) {
			shift = 63 - i*8 -j;
			cout << ((g>>shift)&1);
		}
		cout << endl;
	}
	cout << endl;
	return;
}
