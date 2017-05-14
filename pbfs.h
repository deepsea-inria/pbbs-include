// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "utils.h"
#include "sequence.h"
#include "graph.h"
#include "gettime.h"
namespace pbbs {

using namespace std;

struct nonNegF{bool operator() (intT a) {return (a>=0);}};

pair<intT,intT> pBFS(intT start, graph::graph<intT> GA) {
  timer initTimer;
  initTimer.start();
  intT numVertices = GA.n;
  intT numEdges = GA.m;
  graph::vertex<intT> *G = GA.V;
  intT* Frontier = newA(intT,numEdges);
  intT* Visited = newA(intT,numVertices);
  intT* FrontierNext = newA(intT,numEdges);
  intT* Counts = newA(intT,numVertices);
  {cilk_for(intT i = 0; i < numVertices; i++) Visited[i] = 0;}

  Frontier[0] = start;
  intT frontierSize = 1;
  Visited[start] = 1;

  intT totalVisited = 0;
  int round = 0;
  initTimer.reportStop("initialization");

  timer scanTimer;
  timer mainTimer;
  timer filterTimer;
  while (frontierSize > 0) {
    round++;
    totalVisited += frontierSize;

    scanTimer.start();
    {cilk_for (intT i=0; i < frontierSize; i++) 
      Counts[i] = G[Frontier[i]].degree;}
    intT nr = sequence::scan(Counts,Counts,frontierSize,utils::addF<intT>(),(intT)0);
    scanTimer.stop();

    mainTimer.start();
    // For each vertexB in the frontier try to "hook" unvisited neighbors.
    { cilk_for(intT i = 0; i < frontierSize; i++) {
      intT k= 0;
      intT v = Frontier[i];
      intT o = Counts[i];
      if (G[v].degree < 100) {
        for (intT j=0; j < G[v].degree; j++) {
          intT ngh = G[v].Neighbors[j];
          if (Visited[ngh] == 0 && !__sync_val_compare_and_swap(&Visited[ngh], 0, 1)) {
            FrontierNext[o+j] = ngh;
          } else {
            FrontierNext[o+j] = -1;
          }
        }
      } else {
        cilk_for (intT j=0; j < G[v].degree; j++) {
          intT ngh = G[v].Neighbors[j];
          if (Visited[ngh] == 0 && !__sync_val_compare_and_swap(&Visited[ngh], 0, 1)) {
            FrontierNext[o+j] = ngh;
          } else {
            FrontierNext[o+j] = -1;
          }
        }
      }
     }};
    mainTimer.stop();
    filterTimer.start();    // Filter out the empty solts (marked with -1)
    frontierSize = sequence::filter(FrontierNext,Frontier,nr,nonNegF());
    filterTimer.stop();
  }
  //  std::cerr << totalVisited << " " << round << std::endl;
  scanTimer.reportTotal("scan total");
  mainTimer.reportTotal("main total");
  filterTimer.reportTotal("filter total");
  free(FrontierNext); free(Frontier); free(Counts); free(Visited);
  return pair<intT,intT>(totalVisited,round);
}
} //end namespace
