/*
  File: blurfilter.h

  Declaration of pixel structure and blurfilter function.

 */

#ifndef _BLURFILTER_H_
#define _BLURFILTER_H_
#include <pthread.h>

#ifdef __APPLE__
    #include "barrier.h"
#endif

#include "mpi_data_types.h"

void blurfilter(const int xsize, const int ysize, const int partitioned_height, pixel_t* src, const int radius,
  const double *w, const int my_id, const int n_tasks, pthread_barrier_t* x_done_barrier, pthread_barrier_t* y_done_barrier);

#endif
