#ifndef STOPWATCH_HPP
#define STOPWATCH_HPP
#include <cstdlib>
#include <ctime>

namespace putils {

struct Stopwatch
{
    double acc;
    struct timespec ts;
    struct timespec tf;
    Stopwatch():acc(0.) {}    
    constexpr void clear() noexcept { acc = 0.; }
    void start() { clock_gettime(CLOCK_MONOTONIC,&ts);}
    void stop() { 
        clock_gettime(CLOCK_MONOTONIC,&tf);
        acc += (tf.tv_sec - ts.tv_sec) + 1.e-9 * ( tf.tv_nsec - ts.tv_nsec);
    }
    constexpr double elapsed_time() const noexcept { return acc;}
};

}
#endif