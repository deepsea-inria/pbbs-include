// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2010 Guy Blelloch and Harsha Vardhan Simhadri and the PBBS team
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

#include <iostream>
#include <algorithm>
#include "parallel.h"
#include "utils.h"
#include "sequence.h"
#include "gettime.h"
#include <math.h>
#include "quickSort.h"
#include "transpose.h"

namespace pbbs {
template<class E, class BinPred, class intT>
  void split_positions(E* a, E* b, intT* c, intT length_a, intT length_b, BinPred compare) {
  if (length_a == 0 || length_b == 0) {
    return;
  }
  int pos_a = 0;
  int pos_b = 0;
  int pos_c = 0;
  for (intT i = 0; i < 2 * length_b; i++) {
    c[i] = 0;
  }
  while (pos_b < length_b) {
    while (pos_a < length_a && compare(a[pos_a], b[pos_b])) {
      c[pos_c]++;
      pos_a++;
    }
    pos_c++;
    while (pos_a < length_a && (compare(a[pos_a], b[pos_b]) ^ compare(b[pos_b], a[pos_a]) ^ true)) {
      c[pos_c]++;
      pos_a++;
    }

    pos_b++;
    pos_c++;
    // The pivots are equal
    while (pos_b < length_b && !compare(b[pos_b - 1], b[pos_b])) {
      pos_b++;
      pos_c += 2;
    }
  }
  c[pos_c] = length_a - pos_a;
}

#define SSORT_THR 100000
#define AVG_SEG_SIZE 2
#define PIVOT_QUOT 2

template<class E, class BinPred, class intT>
void sampleSort (E* A, intT n, BinPred f) {
  if (n < SSORT_THR) {
   quickSort(A, n, f);  
//    std::sort(A, A + n, f);
  } else {
    intT sq = (intT)(pow(n,0.5));
    intT rowSize = sq*AVG_SEG_SIZE;
    intT numR = (intT)ceil(((double)n)/((double)rowSize));
    intT numSegs = (sq-1)/PIVOT_QUOT;
    if (numSegs <= 1) {
      std::sort(A, A + n, f);
      return;
    }

    int overSample = 4;
    intT sampleSetSize = numSegs*overSample;
    E* sampleSet = newA(E,sampleSetSize);
    //cout << "n=" << n << " num_segs=" << numSegs << endl;

    // generate samples with oversampling
    cilk_for (intT j=0; j< sampleSetSize; ++j) {
      intT o = utils::hash(j)%n;
      sampleSet[j] = A[o]; 
    }

    // sort the samples
    quickSort(sampleSet, sampleSetSize, f);
//    std::sort(sampleSet, sampleSet + sampleSetSize, f);

    // subselect samples at even stride
    E* pivots = newA(E,numSegs-1);
    cilk_for (intT k=0; k < numSegs-1; ++k) {
      intT o = overSample*k;
      pivots[k] = sampleSet[o];
    }
    free(sampleSet);  
    //nextTime("samples");

    int pivots_size = numSegs - 1;

    numSegs = 2 * numSegs - 1;
    E *B = newA(E, numR*rowSize);
    intT *segSizes = newA(intT, numR*numSegs);
    intT *offsetA = newA(intT, numR*numSegs);
    intT *offsetB = newA(intT, numR*numSegs);

    // sort each row and merge with samples to get counts
    cilk_for (intT r = 0; r < numR; ++r) {
      intT offset = r * rowSize;
      intT size =  (r < numR - 1) ? rowSize : n - offset;
      sampleSort(A+offset, size, f);
      split_positions(A + offset, pivots, segSizes + r*numSegs, size, pivots_size, f);
    }
    //nextTime("sort and merge");

    // transpose from rows to columns
    sequence::scan(segSizes, offsetA, numR*numSegs, plus<intT>(),(intT)0);
    transpose<intT,intT>(segSizes, offsetB).trans(numR, numSegs);
    sequence::scan(offsetB, offsetB, numR*numSegs, plus<intT>(),(intT)0);
    blockTrans<E,intT>(A, B, offsetA, offsetB, segSizes).trans(numR, numSegs);
    {cilk_for (intT i=0; i < n; ++i) A[i] = B[i];}
    //nextTime("transpose");

    free(B); free(offsetA); free(segSizes);

    // sort the columns
    {cilk_for (intT i=0; i<pivots_size + 1; ++i) {
	intT offset = offsetB[(2 * i)*numR];
	if (i == 0) {
	  sampleSort(A, offsetB[numR], f); // first segment
	} else if (i < pivots_size) { // middle segments
	  // if not all equal in the segment
	  if (f(pivots[i-1],pivots[i])) 
	    sampleSort(A+offset, offsetB[(2 * i+1)*numR] - offset, f);
	} else { // last segment
	  sampleSort(A+offset, n - offset, f);
	}
      }
    }
    //nextTime("last sort");
    free(pivots); free(offsetB);
  }
}

#undef compSort
#define compSort(__A, __n, __f) (sampleSort(__A, __n, __f))
} //end namespace
