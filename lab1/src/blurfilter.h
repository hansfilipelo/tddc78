/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.

 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_
#include "blur_mpi_data_types.h"

void blurfilter(const int xsize, const int ysize, pixel_t* src, const int radius,
  const double *w, const int my_rank, const int n_tasks);

#endif
