#include "thresfilter.h"
#define uint unsigned int

double calculate_partitioned_mean(const int nr_pixels, const unsigned total_pixels, const pixel_t* src, const unsigned my_id)
{
    uint i_start, i_stop;
    unsigned int sum = 0;
    unsigned int i;
    double npix = (double)nr_pixels;

    if (my_id != _N_TASKS_) {
        i_start = my_id*nr_pixels;
        i_stop = my_id*nr_pixels+nr_pixels;
    }
    else{
        i_start = my_id*nr_pixels+nr_pixels;
        i_stop = total_pixels;
    }

    for(i = i_start, sum = 0; i < i_stop; i++) {
        sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    }

    if (nr_pixels != 0) {
        return sum / npix;
    }

    return 0;
}

// ---------

unsigned int calculate_total_mean(double* means, const unsigned int nr_of_means, const unsigned int partitioned_pixels, const unsigned int remainder_pixels, const unsigned int total_pixels)
{
    double total_mean = 0;
    unsigned int i;

    for (i = 0; i < nr_of_means-1; i++) {
        total_mean += means[i]*partitioned_pixels;
    }
    total_mean += means[i]*remainder_pixels;

    total_mean /= total_pixels;

    return (uint)total_mean;
}

// ---------

void thresfilter(const int nr_pixels, const unsigned total_pixels, const unsigned int total_mean, pixel_t* src, const unsigned my_id){

    uint pixel_sum, i, i_start, i_stop;

    if (my_id != _N_TASKS_) {
        i_start = my_id*nr_pixels;
        i_stop = my_id*nr_pixels+nr_pixels;
    }
    else{
        i_start = my_id*nr_pixels+nr_pixels;
        i_stop = total_pixels;
    }

    for(i = i_start ; i < i_stop; i++) {
        pixel_sum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
        if(total_mean > pixel_sum) {
            src[i].r = src[i].g = src[i].b = 0;
        }
        else {
            src[i].r = src[i].g = src[i].b = 255;
        }
    }
}
