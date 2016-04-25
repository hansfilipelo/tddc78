/*
  File: thresfilter.h

  Declaration of pixel structure and thresfilter function.

 */
#ifndef _THRESFILTER_H_
#define _THRESFILTER_H_

#include "mpi_data_types.h"

double calculate_partitioned_mean(const int nr_pixels, const pixel_t* src);
unsigned int calculate_total_mean(double* means, const unsigned int nr_of_means, const unsigned int partitioned_pixels, const unsigned int remainder_pixels, const unsigned int total_pixels);
void thresfilter(const int nr_pixels, const unsigned int total_mean, pixel_t* src);

#endif
