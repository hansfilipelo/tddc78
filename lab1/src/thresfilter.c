#include "thresfilter.h"
#define uint unsigned int

double calculate_partitioned_mean(const int nr_pixels, const pixel_t* src)
{
  unsigned int sum, i;
  double npix = (double)nr_pixels;

  for(i = 0, sum = 0; i < nr_pixels; i++) {
    sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
  }

  return sum / npix;
}

// ---------

unsigned int calculate_total_mean(double* means, const unsigned int nr_of_means, const unsigned int partitioned_pixels, const unsigned int remainder_pixels, const unsigned int total_pixels)
{
  double total_mean = 0;
  unsigned int i;

  for (i = 0; i < nr_of_means-1; i++) {
    total_mean += means[i]*partitioned_pixels/total_pixels;
  }
  total_mean += means[i]*remainder_pixels/total_pixels;

  return (uint)total_mean;
}

// ---------

void thresfilter(const int nr_pixels, const unsigned int total_mean, pixel_t* src){

  uint pixel_sum, i;

  for(i = 0; i < nr_pixels; i++) {
    pixel_sum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(total_mean > pixel_sum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }
}
