#include <mpi.h>
#ifndef _MPI_DATA_TYPES_H_
#define _MPI_DATA_TYPES_H_

typedef struct {
  unsigned int width;
  unsigned int height;
} size_data_t;

typedef struct {
  unsigned char r;
  unsigned char g;
  unsigned char b;
} pixel_t;

void create_mpi_size_data(const size_data_t* size_data, MPI_Datatype* mpi_size_data);

void create_mpi_pixel(const pixel_t* pixel, MPI_Datatype* mpi_pixel);

#endif // _MPI_DATA_TYPES_H_
