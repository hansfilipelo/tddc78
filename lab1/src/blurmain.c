#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include "blur_mpi_data_types.h"

#define MAX_RAD 1000

int main (int argc, char ** argv) {
  int n_tasks, my_rank;
  int radius, xsize, ysize, colmax;
  double w[MAX_RAD];
  pixel_t* src;

  /* MPI initilization */
  MPI_Comm com = MPI_COMM_WORLD;
  MPI_Init( &argc, &argv );
  MPI_Comm_size( com, &n_tasks );
  MPI_Comm_rank( com, &my_rank );

  /* Take care of the arguments */
  if (argc != 4) {
    if(my_rank == 0) {
      fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
    }
    exit(1);
  }
  radius = atoi(argv[1]);
  if((radius > MAX_RAD) || (radius < 1)) {
    if(my_rank == 0) {
      fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
    }
    exit(1);
  }

  /* If my rank is zero read input image */
  if(my_rank == 0)
  {
    src = malloc(sizeof(pixel_t) * MAX_PIXELS);
    struct timespec stime, etime; // TODO: Mac

    /* read file */
    if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
      exit(1); // TODO: Broadcast exit

    if (colmax > 255) {
      fprintf(stderr, "Too large maximum color-component value\n");
      exit(1); // TODO: Broadcast exit
    }

    printf("Has read the image, generating coefficients\n");
  }

  // TODO: Create types

  // TODO: Broadcast width and heigth

  // TODO: Scatter pixels

  /* Calculate gaussian weights */
  get_gauss_weights(radius, w);

  printf("Calling filter\n");

  #ifdef _linux_
  clock_gettime(CLOCK_REALTIME, &stime); // TODO: Mac
  #endif

  blurfilter(xsize, ysize, src, radius, w); // TODO: Just my part

  #ifdef _linux_
  clock_gettime(CLOCK_REALTIME, &etime); // TODO: Mac

  // TODO: Gather pixels

  // TODO: Time only for rank 0 from scatter to gather
  printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
  1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;
  #endif

  /* write result */
  printf("Writing output file\n");

  // TODO: Just rank 0
  if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0)
  exit(1);


  return(0);
}
