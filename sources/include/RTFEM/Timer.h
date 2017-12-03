#ifndef PROJECT_TIMER_H
#define PROJECT_TIMER_H

#include <chrono>

namespace rtfem {

class Timer {
public:
    void Start() {
        start_time_ = std::chrono::system_clock::now();
    }

    double Stop() {
        auto end_time = std::chrono::system_clock::now();
        elapsed_ = end_time - start_time_;

        elapsed_computed_ = true;

        return elapsed_.count() * 1000;
    }

    double Elapsed() const {
        if(!elapsed_computed_){
            return 0;
        }

        return elapsed_.count();
    }

    double ElapsedMS() const {
        if(!elapsed_computed_)
            return 0;
        return elapsed_.count() * 1000;
    }
private:
    std::chrono::system_clock::time_point start_time_;
    std::chrono::duration<double> elapsed_;

    bool elapsed_computed_ = false;
};
}


#endif //PROJECT_TIMER_H
