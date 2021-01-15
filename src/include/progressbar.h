#ifndef PROGRESS_BAR_H
#define PROGRESS_BAR_H

#include <string>
#include <sys/time.h>

using namespace std;

class ProgressBar {



private:
    int size;
    double start_time;
    string description;
    float prv_val;

    struct timeval start;
    string make_bar(float);
    double get_time();

public:

    ProgressBar(int);

    void update (float, string);
};

#endif //PROGRESS_BAR_H
