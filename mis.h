#include "sequence.h"
#include "graph.h"
#include "parallel.h"
#include "speculativefor.h"

namespace pbbs {

using namespace std;

struct MISstep {
  char flag;
  char *Flags;  graph::vertex<intT>*G;
  MISstep(char* _F, graph::vertex<intT>* _G) : Flags(_F), G(_G) {}

  bool reserve(intT i) {
    intT d = G[i].degree;
    flag = 1;
    for (intT j = 0; j < d; j++) {
      intT ngh = G[i].Neighbors[j];
      if (ngh < i) {
	if (Flags[ngh] == 1) { flag = 2; return 1;}
	// need to wait for higher priority neighbor to decide
	else if (Flags[ngh] == 0) flag = 0; 
      }
    }
    return 1;
  }

  bool commit(intT i) { return (Flags[i] = flag) > 0;}
};

char* maximalIndependentSet(graph::graph<intT> GS) {
  intT n = GS.n;
  graph::vertex<intT>* G = GS.V;
  char* Flags = newArray(n, (char) 0);
  MISstep mis(Flags, G);
  speculative_for(mis, 0, n, 20);
  return Flags;
}

} // end namespace