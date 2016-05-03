#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "ppmio.h"
#include "thresfilter.h"
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
  int n_tasks, my_rank, colmax, total_pixels;

  struct timespec stime, etime;

  size_data_t size_data;
  unsigned int partitioned_pixels;
  unsigned int remainder_pixels;
  pixel_t* src;
  double* all_partitioned_means;
  unsigned int total_mean;

  /* MPI initilization */
  MPI_Comm com = MPI_COMM_WORLD;
  MPI_Init( &argc, &argv );
  MPI_Comm_size( com, &n_tasks );
  MPI_Comm_rank( com, &my_rank );

  /* Take care of the arguments */
  if (argc != 3) {
    if(my_rank == 0) {
      fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
    }
    exit(1);
  }

  /* If my rank is zero read input image */
  if(my_rank == 0)
  {
    src = malloc(sizeof(pixel_t) * MAX_PIXELS);
    all_partitioned_means = malloc(sizeof(double) * (n_tasks+1));

    /* read file */
    if(read_ppm (argv[1], &size_data.width, &size_data.height, &colmax, (char *) src) != 0) {
      exit(1); // TODO: Broadcast exit
    }

    if (colmax > 255) {
      fprintf(stderr, "Too large maximum color-component value\n");
      exit(1); // TODO: Broadcast exit
    }

    printf("Has read the image.\n");
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
  total_pixels = size_data.height*size_data.width;
  partitioned_pixels = total_pixels/n_tasks;
  remainder_pixels = total_pixels-(n_tasks*partitioned_pixels);
  double partitioned_mean;

  if(my_rank != 0) {
      src = malloc(sizeof(pixel_t) * partitioned_pixels);
  }

  if(my_rank == 0) {
    printf("Scatter picture to filter it\n");
    clock_gettime(CLOCK_REALTIME, &stime);
  }

  // Scatter data to different processes
  if (my_rank == 0) {
    MPI_Scatter(src, partitioned_pixels, mpi_pixel, MPI_IN_PLACE, partitioned_pixels, mpi_pixel, 0, com);
  }
  else{
    MPI_Scatter(src, partitioned_pixels, mpi_pixel, src, partitioned_pixels, mpi_pixel, 0, com);
  }

  partitioned_mean = calculate_partitioned_mean(partitioned_pixels, src);

  MPI_Gather(&partitioned_mean, 1, MPI_DOUBLE, all_partitioned_means, 1, MPI_DOUBLE, 0, com);

  if(my_rank == 0)
  {
    //all_partitioned_means[0] = partitioned_mean; // Not needed since gather will put it here unless MPI_IN_PLACE specified
    all_partitioned_means[n_tasks] = calculate_partitioned_mean(remainder_pixels, &src[total_pixels - remainder_pixels]);
    total_mean = calculate_total_mean(all_partitioned_means, n_tasks+1, partitioned_pixels, remainder_pixels, total_pixels);
  }

  MPI_Bcast(&total_mean, 1, MPI_DOUBLE, 0, com);

  thresfilter(partitioned_pixels, total_mean, src);

  if (my_rank == 0) {
    thresfilter(remainder_pixels, total_mean, &src[total_pixels-remainder_pixels]);
  }

  if (my_rank == 0) {
    MPI_Gather(MPI_IN_PLACE, partitioned_pixels, mpi_pixel, src, partitioned_pixels, mpi_pixel, 0, com);
  }
  else{
    MPI_Gather(src, partitioned_pixels, mpi_pixel, src, partitioned_pixels, mpi_pixel, 0, com);
  }

  if(my_rank == 0) {
    clock_gettime(CLOCK_REALTIME, &etime);
    printf("Gathered picture after filtering\n");

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
    1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");

    if(write_ppm (argv[2], size_data.width, size_data.height, (char *)src) != 0){
      exit(1);
    }
  }

  MPI_Finalize();
  return(0);
}
