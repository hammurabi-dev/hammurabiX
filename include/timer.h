/** 
 * timer module
 * timer record time elapse from Timer::start to Timer::stop
 * timer support nested/overlaped recording
 */
#ifndef HAMMURABI_TIMR_H
#define HAMMURABI_TIMR_H

#include <iostream>
#include <chrono>
#include <ctime>
#include <string>
#include <map>
#include <cassert>

class Timer{
    typedef std::chrono::time_point<std::chrono::high_resolution_clock> tick;
    typedef std::chrono::duration<double,std::milli> duration;
    typedef std::map<std::string,std::pair<tick,double>> timecache;
#ifndef NDEBUG
public:
#endif
    timecache record;
public:
    /**
     * record start point of timer
     * 1st argument: name of timing record
     */
    inline void start(std::string flag){
        record[flag].first = std::chrono::high_resolution_clock::now();
    }
    /**
     * record end point of timer
     * then calculate corresponding elapsed time
     * 1st argument: name of timing record
     */
    inline void stop(std::string flag){
        duration diff = (std::chrono::high_resolution_clock::now() - record[flag].first);
        record[flag].second = diff.count();
    }
    /**
     * print to stdout the elapsed time in ms resolution
     * 1st argument: (optional) name of timing record
     * if no name provided, print all records
     */
    inline void print(std::string flag=std::string()){
        if(flag.empty()){
            for(auto i: record)
                std::cout<<i.first<<" "<<i.second.second<<" ms"<<std::endl;
        }
        else{
            auto target = record.find(flag);
            assert(target!=record.end());
            std::cout<<target->first<<" "<<target->second.second<<" ms"<<std::endl;
        }
    }
};

#endif
