#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"
#include "mpi_data_types.h"

#ifdef __APPLE__
  #include "timing_mach.h"
#endif

#ifdef _linux_
  #include <sys/time.h>
  #include <time.h>
#endif

int main (int argc, char ** argv) {
  int n_tasks, my_rank;
  int radius, colmax;
  int offset;

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

  /* Check that nr of tasks is larger than one */
  if(n_tasks < 2) {
    fprintf(stderr, "This implementation doesn't support fewer than two tasks\n");
    exit(1);
  }

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

  double* w = malloc(sizeof(double)*radius);

  /* If my rank is zero read input image */
  if(my_rank == 0)
  {
    src = malloc(sizeof(pixel_t) * MAX_PIXELS);

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
  int send_count[n_tasks];
  int receive_count[n_tasks];
  int displacements[n_tasks];
  int receive_displacements[n_tasks];

  // Since the ranks needs different amounts of data, we need to specify for
  // MPI_Scatterv (<-- "v") which addresses each rank need to read from
  send_count[0] = (partitioned_height + radius)*size_data.width;
  displacements[0] = 0;
  receive_count[0] = partitioned_height*size_data.width;
  receive_displacements[0] = 0;

  int i = 0;
  for (i = 1; i < n_tasks-1; i++) {
    send_count[i] = (2*radius + partitioned_height)*size_data.width;
    receive_count[i] = partitioned_height*size_data.width;
    displacements[i] = (partitioned_height*i-radius)*size_data.width - 1;
    receive_displacements[i] = (partitioned_height*i)*size_data.width - 1;
  }
  send_count[n_tasks-1] = (radius + partitioned_height + remainder_height)*size_data.width;
  receive_count[n_tasks-1] = (partitioned_height+remainder_height)*size_data.width;
  displacements[n_tasks-1] = (partitioned_height*(n_tasks-1)-radius)*size_data.width - 1;
  receive_displacements[n_tasks-1] = (partitioned_height*(n_tasks-1))*size_data.width - 1;


  if(my_rank != 0) {
      src = malloc(sizeof(pixel_t) * send_count[my_rank]);
  }
  // // TEST TEST TEST
  // if(my_rank == 0) {
  //   printf("partitioned_height = %i, remainder_height = %i \n", partitioned_height, remainder_height);
  //   printf("Width = %i, height = %i, send_count(my_rank) = %i, displacements(my_rank) = %i \n", size_data.width, size_data.height, send_count[my_rank], displacements[my_rank]);
  // }

  /* Calculate gaussian weights */
  get_gauss_weights(radius, w);

  if(my_rank == 0) {
    printf("Scatter picture to filter it\n");
    clock_gettime(CLOCK_REALTIME, &stime);
  }

  // Scatter data to different processes
  if (my_rank == 0) {
    MPI_Scatterv(src, send_count, displacements, mpi_pixel, MPI_IN_PLACE, send_count[my_rank], mpi_pixel, 0, com);
  }
  else{
    MPI_Scatterv(src, send_count, displacements, mpi_pixel, src, send_count[my_rank], mpi_pixel, 0, com);
  }

  // Run the filter
  blurfilter(size_data.width, send_count[my_rank]/size_data.width, src, radius, w, my_rank, n_tasks);

  // Gather data to rank 0
  if (my_rank == 0) {
    MPI_Gatherv(MPI_IN_PLACE, receive_count[my_rank], mpi_pixel, src, receive_count, receive_displacements, mpi_pixel, 0, com);
  }
  else {
    MPI_Gatherv(&src[radius*size_data.width], receive_count[my_rank], mpi_pixel, src, receive_count, receive_displacements, mpi_pixel, 0, com);
  }

  if(my_rank == 0) {
    clock_gettime(CLOCK_REALTIME, &etime);
    printf("Gather picture after filtering\n");

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
    1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");

    if(write_ppm (argv[3], size_data.width, size_data.height, (char *)src) != 0){
      exit(1);
    }
  }

  MPI_Finalize();
  return(0);
}
