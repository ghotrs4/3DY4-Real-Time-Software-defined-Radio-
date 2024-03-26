#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"

#include <iostream>
#include <thread>
#include <cstdio>
#include <queue>
#include <mutex>
#include <condition_variable>

using namespace std;


class threadSafeQ{

public:

    std::queue<std::vector<float>> q;
    std::mutex m;
    std::condition_variable c;

  //Adds element to queue (with thread safe implementation)
  void enqueue(std::vector<float> v) {
    std::unique_lock<std::mutex> lock(m);
    q.push(v);
    c.notify_one();
  }

   //Gets element from queue (with thread safe implementation)
    std::vector<float> dequeue(void) {
    std::unique_lock<std::mutex> lock(m);
    
    //Waits until queue is not empty
    while(q.empty()) {
      //Releases lock when notified by enqueue
      c.wait(lock);
    }
    
    //Defines the value to be returned as the front element of the queue
    std::vector<float> val = q.front();
    q.pop();
    return val;
  }
  
  //Checks if queue is empty
  bool isEmpty(void) {
    return this->q.empty();
      
  }

};
