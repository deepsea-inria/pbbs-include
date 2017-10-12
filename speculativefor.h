#include "parallel.h"
#include "utils.h"
#include "sequence.h"

namespace pbbs {

struct reservation {
  intT r;
  reservation() : r(INT_T_MAX) {}
  void reserve(intT i) { utils::writeMin(&r, i); }
  bool reserved() { return (r < INT_T_MAX);}
  void reset() {r = INT_T_MAX;}
  bool check(intT i) { return (r == i);}
  bool checkReset(intT i) {
    if (r==i) { r = INT_T_MAX; return 1;}
    else return 0;
  }
};

inline void reserveLoc(intT& x, intT i) {utils::writeMin(&x,i);}

template <class S>
intT speculative_for(S step, intT s, intT e, int granularity, 
		     bool hasState = 1, int maxTries = -1) {
  if (maxTries < 0) maxTries = 100 + 200 * granularity;
  intT maxRoundSize = (e - s) / granularity + 1;
  intT* I = newA(intT,maxRoundSize);
  intT* Ihold = newA(intT, maxRoundSize);
  bool* keep = newA(bool, maxRoundSize);
  S *state;
  if (hasState) {
    state = newA(S, maxRoundSize);
    for (intT i = 0; i < maxRoundSize; i++) state[i] = step;
  }

  int round = 0; 
  intT numberDone = s; // number of iterations done
  intT numberKeep = 0; // number of iterations to carry to next round
  intT totalProcessed = 0;

  while (numberDone < e) {
    //cout << "numberDone=" << numberDone << endl;
    if (round++ > maxTries) {
      std::cout << "speculativeLoop: too many iterations, increase maxTries parameter" << std::endl;
      abort();
    }
    intT size = std::min(maxRoundSize, e - numberDone);
    totalProcessed += size;

    if (hasState) {
      cilk_for (intT i = 0; i < size; i++) {
	if (i >= numberKeep) I[i] = numberDone + i;
	keep[i] = state[i].reserve(I[i]);
      } 
    } else {
      cilk_for (intT i = 0; i < size; i++) {
	if (i >= numberKeep) I[i] = numberDone + i;
	keep[i] = step.reserve(I[i]);
      } 
    }

    if (hasState) {
      cilk_for (intT i = 0; i < size; i++) 
	if (keep[i]) keep[i] = !state[i].commit(I[i]);
    } else {
      cilk_for (intT i = 0; i < size; i++) 
	if (keep[i]) keep[i] = !step.commit(I[i]);
    }

    // keep edges that failed to hook for next round
    numberKeep = sequence::pack(I, Ihold, keep, size);
    std:: swap(I, Ihold);
    numberDone += size - numberKeep;
  }
  free(I); free(Ihold); free(keep); 
  if(hasState) free(state);
  return totalProcessed;
}


} // end namespace
