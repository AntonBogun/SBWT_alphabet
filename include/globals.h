#pragma once
#include <chrono>
#define __get_now() std::chrono::high_resolution_clock::now()
#define __get_duration(text, start, end) std::cout << text << " took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl; 
// using namespace std;


class ProgressPrinter {
public:
    int64_t n_jobs;
    int64_t processed;
    int64_t total_prints;
    int64_t next_print;
    bool first_print;
    std::chrono::time_point<std::chrono::high_resolution_clock> start;

    ProgressPrinter(int64_t n_jobs, int64_t total_prints) 
        : n_jobs(n_jobs), processed(0), total_prints(total_prints), next_print(0), first_print(true) {}

    void jobDone() {
        if(first_print){
            start = __get_now();
        }
        if(next_print == processed) {
            if(!first_print) std::cerr << '\r' << std::flush; // Erase current line
            first_print = false;
            
            int64_t progress_percent = std::round(100 * ((double)processed / n_jobs));
            std::cerr << std::to_string(progress_percent) + "%, " ;
            // expected seconds left
            auto now = __get_now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
            auto expected_duration = duration * n_jobs / processed;
            auto expected_seconds_left = (expected_duration - duration) / 1000;
            std::cerr << "ETA: " << expected_seconds_left << "s" << std::flush;
            next_print += n_jobs / (4*total_prints);
        }
        processed++;
        if(processed == n_jobs) std::cerr << "\r100%" << std::endl; // No more prints coming
    }
};