#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include "blur_mpi_data_types.h"

#ifdef __APPLE__
#include "timing_mach.h"
#endif

#define MAX_RAD 1000

int main (int argc, char ** argv) {
  int n_tasks, my_rank;
  int radius, colmax;
  double w[MAX_RAD];

  struct timespec stime, etime;

  size_data_t size_data;
  unsigned int partitioned_height;
  unsigned int remainder_height;
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
    if(read_ppm (argv[2], &size_data.width, &size_data.height, &colmax, (char *) src) != 0) {
      exit(1); // TODO: Broadcast exit
    }

    if (colmax > 255) {
      fprintf(stderr, "Too large maximum color-component value\n");
      exit(1); // TODO: Broadcast exit
    }

    printf("Has read the image, generating coefficients\n");
  }

  // Create data types for MPI
  MPI_Datatype mpi_size_data;
  MPI_Datatype mpi_pixel;
  pixel_t dummy_pixel;
  create_mpi_size_data(&size_data, &mpi_size_data);
  create_mpi_pixel(&dummy_pixel, &mpi_pixel);

  // Broadcast width and heigth
  MPI_Bcast(&size_data, 1, mpi_size_data, 0, com);

  // Calculate some sweet height values
  partitioned_height = size_data.height/n_tasks;
  remainder_height = size_data.height-(n_tasks*partitioned_height);

  // Now we need to partition the task correctly for MPI_scatterv
  int send_counts[n_tasks];
  int displacements[n_tasks];

  // Since the ranks needs different amounts of data, we need to specify for
  // MPI_Scatterv (<-- "v") which addresses each rank need to read from
  displacements[0] = 0;
  send_counts[0] = partitioned_height+radius;

  for (size_t i = 1; i < n_tasks-1; i++) {
    send_counts[i] = 2*radius+partitioned_height;
    displacements[i] = partitioned_height*i-radius-1;
  }
  send_counts[n_tasks-1] = radius+partitioned_height+remainder_height;
  displacements[n_tasks-1] = partitioned_height*(n_tasks-1)-radius-1;

  // Scatter data to different processes
  MPI_Scatterv(src, send_counts, displacements, mpi_pixel, src, send_counts[my_rank], mpi_pixel, 0, com);

  /* Calculate gaussian weights */
  get_gauss_weights(radius, w);

  printf("Calling filter\n");

  clock_gettime(CLOCK_REALTIME, &stime); // TODO: Mac

  blurfilter(size_data.width, send_counts[my_rank], src, radius, w); // TODO: Just my part

  MPI_Gatherv(src, send_counts[my_rank], mpi_pixel, src, send_counts, displacements, mpi_pixel, 0, com);

  clock_gettime(CLOCK_REALTIME, &etime); // TODO: Mac


  // TODO: Time only for rank 0 from scatter to gather
  printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
  1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

  /* write result */
  printf("Writing output file\n");

  // TODO: Just rank 0
  if (my_rank == 0){
    if(write_ppm (argv[3], size_data.width, size_data.height, (char *)src) != 0){
      exit(1);
    }
  }

  return(0);
}
