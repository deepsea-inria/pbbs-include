// -*- C++ -*-

#ifndef _BENCH_GETTIME_INCLUDED
#define _BENCH_GETTIME_INCLUDED

#include <stdlib.h>
#include <sys/time.h>
#include <iomanip>
#include <iostream>

namespace pbbs {

struct timer {
  double totalTime;
  double lastTime;
  double totalWeight;
  bool on;
  struct timezone tzp;
  timer() {
    struct timezone tz = {0, 0};
    totalTime=0.0; 
    totalWeight=0.0;
    on=0; tzp = tz;}
  double getTime() {
#ifdef TIME_MEASURE
    timeval now;
    gettimeofday(&now, &tzp);
    return ((double) now.tv_sec) + ((double) now.tv_usec)/1000000.;
#endif
    return 0.0;
  }
  void start () {
#ifdef TIME_MEASURE
    on = 1;
    lastTime = getTime();
#endif
  } 
  double stop () {
#ifdef TIME_MEASURE
    on = 0;
    double d = (getTime()-lastTime);
    totalTime += d;
    return d;
#endif
    return 0.0;
  }
  double stop (double weight) {
#ifdef TIME_MEASURE
    on = 0;
    totalWeight += weight;
    double d = (getTime()-lastTime);
    totalTime += weight*d;
    return d;
#endif
    return 0.0;
  }

  double total() {
#ifdef TIME_MEASURE
    if (on) return totalTime + getTime() - lastTime;
    else return totalTime;
#endif
    return 0.0;
  }

  double next() {
#ifdef TIME_MEASURE
    if (!on) return 0.0;
    double t = getTime();
    double td = t - lastTime;
    totalTime += td;
    lastTime = t;
    return td;
#endif
    return 0.0;
  }

  void reportT(double time) {
#ifdef TIME_MEASURE
    std::cout << "PBBS-time: " << std::setprecision(3) << time <<  std::endl;;
#endif
  }

  void reportTime(double time) {
#ifdef TIME_MEASURE
    reportT(time);
#endif
  }

  void reportStop(double weight, std::string str) {
#ifdef TIME_MEASURE
    std::cout << str << " :" << weight << ": ";
    reportTime(stop(weight));
#endif
  }

  void reportStop(std::string str) {
#ifdef TIME_MEASURE
    std::cout << str << ": ";
    reportTime(stop());
#endif
  }

  void reportTotal() {
#ifdef TIME_MEASURE
    double to = (totalWeight > 0.0) ? total()/totalWeight : total();
    reportTime(to);
    totalTime = 0.0;
    totalWeight = 0.0;
#endif
  }

  void reportTotal(std::string str) {
#ifdef TIME_MEASURE
    std::cout << str << " : "; 
    reportTotal();
#endif
  }

  void reportNext() {
#ifdef TIME_MEASURE
    reportTime(next());
#endif
  }

  void reportNext(std::string str) {std::cout << str << " : "; reportNext();}
};

static timer _tm;
#ifdef TIME_MEASURE
#define timeStatement(_A,_string) _tm.start();  _A; _tm.reportNext(_string);
#define startTime() _tm.start();
#define stopTime(_weight,_str) _tm.reportStop(_weight,_str);
#define reportTime(_str) _tm.reportTotal(_str);
#define nextTime(_string) _tm.reportNext(_string);
#define nextTimeN() _tm.reportT(_tm.next());
#else
#define timeStatement(_A,_string)
#define startTime() 
#define stopTime(_weight,_str) 
#define reportTime(_str) 
#define nextTime(_string) 
#define nextTimeN() _tm.reportT(_tm.next())
#endif

} //end namespace
#endif // _BENCH_GETTIME_INCLUDED

