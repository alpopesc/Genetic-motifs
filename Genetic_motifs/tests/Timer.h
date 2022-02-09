

#ifndef TEAM8_TIMER_H
#define TEAM8_TIMER_H

#include <chrono>
#include <iostream>
#include <thread>
#include <memory>

/*! \class Network
  A Timer measures the the time of the scope in which it was initialised and prints the time on the terminal.

 */

class Timer{

public:
    Timer(){
        m_StartTimepoint = std::chrono::high_resolution_clock::now();
    }

    ~Timer(){
        Stop();
    }

    void Stop() {
        auto endTimepoint = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch();

        auto duration = end - start;
        double ms = duration.count()*0.000001;

        std::cout << ms << "s";
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock > m_StartTimepoint;

};


#endif 
