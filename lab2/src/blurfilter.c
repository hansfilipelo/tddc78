/*
File: blurfilter.c

Implementation of blurfilter function.

*/
#include <stdio.h>
#include <stdlib.h>
#include "blurfilter.h"
#include "ppmio.h"


pixel_t* pix(pixel_t* image, const int xx, const int yy, const int xsize)
{
  register int off = xsize*yy + xx;

  #ifdef DBG
  if(off >= MAX_PIXELS) {
    fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
  }
  #endif
  return (image + off);
}

void blurfilter(const int xsize, const int ysize, const int partitioned_height, pixel_t* src, const int radius,
  const double *w, const int my_id, const int n_tasks, pthread_barrier_t* x_done_barrier, pthread_barrier_t* y_done_barrier)
  {
    int x,y,x2,y2, wi;
    double r,g,b,n, wc;
    pixel_t* dst = malloc(sizeof(pixel_t)*xsize*ysize);

    // Set start and stop values
    int y_start = my_id*partitioned_height, y_stop = my_id*partitioned_height+partitioned_height;
    int x_start = 0, x_stop = xsize;

    if (my_id != 0) {
        y_start -= radius;
    }
    if (my_id != n_tasks-1) {
        y_stop += radius;
    }
    if ( my_id == n_tasks-1 ) {
        y_stop = ysize;
    }

    for (y=y_start; y<y_stop; y++) {
      for (x=x_start; x<x_stop; x++) {
        r = w[0] * pix(src, x, y, xsize)->r;
        g = w[0] * pix(src, x, y, xsize)->g;
        b = w[0] * pix(src, x, y, xsize)->b;
        n = w[0];
        for ( wi=1; wi <= radius; wi++) {
          wc = w[wi];
          x2 = x - wi;
          if(x2 >= 0) {
            r += wc * pix(src, x2, y, xsize)->r;
            g += wc * pix(src, x2, y, xsize)->g;
            b += wc * pix(src, x2, y, xsize)->b;
            n += wc;
          }
          x2 = x + wi;
          if(x2 < xsize) {
            r += wc * pix(src, x2, y, xsize)->r;
            g += wc * pix(src, x2, y, xsize)->g;
            b += wc * pix(src, x2, y, xsize)->b;
            n += wc;
          }
        }
        pix(dst,x,y, xsize)->r = r/n;
        pix(dst,x,y, xsize)->g = g/n;
        pix(dst,x,y, xsize)->b = b/n;
      }
    }

    // Set different start and stop values depending on current rank
    if(my_id != 0) {
      y_start += radius;
    }
    if(my_id != n_tasks-1) {
      y_stop -= radius;
    }

    pthread_barrier_wait(x_done_barrier);

    for (y=y_start; y<y_stop; y++) {
      for (x=x_start; x<x_stop; x++) {
        r = w[0] * pix(dst, x, y, xsize)->r;
        g = w[0] * pix(dst, x, y, xsize)->g;
        b = w[0] * pix(dst, x, y, xsize)->b;
        n = w[0];
        for ( wi=1; wi <= radius; wi++) {
          wc = w[wi];
          y2 = y - wi;
          if(y2 >= 0) {
            r += wc * pix(dst, x, y2, xsize)->r;
            g += wc * pix(dst, x, y2, xsize)->g;
            b += wc * pix(dst, x, y2, xsize)->b;
            n += wc;
          }
          y2 = y + wi;
          if(y2 < ysize) {
            r += wc * pix(dst, x, y2, xsize)->r;
            g += wc * pix(dst, x, y2, xsize)->g;
            b += wc * pix(dst, x, y2, xsize)->b;
            n += wc;
          }
        }
        pix(src,x,y, xsize)->r = r/n;
        pix(src,x,y, xsize)->g = g/n;
        pix(src,x,y, xsize)->b = b/n;
      }
    }

    pthread_barrier_wait(y_done_barrier);
  }