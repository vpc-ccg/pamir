#include <cstdio>
#include <cmath>

#include "progressbar.h"

ProgressBar::ProgressBar(int size) : size(size) {
    start_time = get_time ();
}

string ProgressBar::make_bar(float progress) {
    int steps = (int)(progress * size / 100.0);
    string out = "[";
    for (int i=0;i<steps; i++)
        out += '=';

    if (progress<100)
        out += ">";

    for (int i=steps+1; i<size; i++)
        out += " ";
    out += "]";

    return out;
}

void ProgressBar::update(float progress, string prefix) {
    if (progress <= 0.0000000000000001 || progress >= 100.0 - 0.0000000000000001 || progress - prv_val > 0.1) {
        string out = make_bar(progress);
        if (progress == 0)
            start_time = get_time();
        double elapsed_time = get_time() - start_time;
        int hh = (int)elapsed_time / 3600;
        int mm = (((int)elapsed_time) % 3600) / 60;
        int ss = (int) elapsed_time % 60;

        fprintf (stderr, "\r%s %s %6.2f%% Elapsed Time: %0.2d:%0.2d:%0.2d", prefix.c_str(), out.c_str(), progress, hh, mm, ss);
        if (abs(100 -progress )< 0.0000000000000001) {
            fprintf(stderr, "\n");
        }
        prv_val = progress;
    }
}

double ProgressBar::get_time(void)
{
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec+t.tv_usec/1000000.0;
}