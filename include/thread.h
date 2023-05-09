#ifndef DY4_THREAD_H
#define DY4_THREAD_H

#include <cmath>
#include <functional>
#include <thread>
#include <condition_variable>
#include <queue>
#include <mutex>

void front_end_producer(std::queue<std::vector<float>> &, std::mutex &, std::condition_variable &, bool &, int );
void stereo_consumer(std::queue<std::vector<float>> &, std::mutex &, std::condition_variable &, bool &,int );
void data(pll_state &samples);

#endif
