#ifndef _PARALLEL_H_
#define _PARALLEL_H_

// intel cilk+
#if defined(USE_CILK_PLUS_RUNTIME) || defined(USE_PASL_RUNTIME)
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
//#define parallel_for cilk_for
#define parallel_main main
#define parallel_for_1 _Pragma("cilk grainsize = 1") cilk_for
#define parallel_for_256 _Pragma("cilk grainsize = 256") cilk_for

// c++
#else
#define cilk_spawn
#define cilk_sync
#define parallel_main main
#define parallel_for for
#define parallel_for_1 for
#define parallel_for_256 for
#define cilk_for for

#endif

#include <limits.h>

#if defined(LONG)
typedef long intT;
typedef unsigned long uintT;
#define INT_T_MAX LONG_MAX
#else
typedef int intT;
typedef int uintT;
#define INT_T_MAX INT_MAX
#endif
#endif