#include "blur_mpi_data_types.h"

// Create mpi size data type
void create_mpi_size_data(size_data_t* size_data, MPI_Datatype* mpi_size_data)
{
  int block_lengths [] = { 1, 1 }; // Lengths of type elements
  MPI_Datatype block_types [] =
  {
    MPI_UNSIGNED_INT,
    MPI_UNSIGNED_INT
  }; //Set types
  MPI_Aint start, displ[2];
  
  MPI_Address( &size_data, &start );
  MPI_Address( &size_data->width, &displ[0] );
  MPI_Address( &size_data->height, &displ[1] );

  displ[0] -= start; // Displacement relative to address of start
  displ[1] -= start; // Displacement relative to address of start

  MPI_Type_struct( 2, block_lengths, displ, block_types, &mpi_size_data );
  MPI_Type_commit( &mpi_size_data );
}

// Create mpi pixel data type
void create_mpi_pixel(pixel_t* pixel, MPI_Datatype* mpi_pixel)
{
  int block_lengths [] = { 1, 1, 1 }; // Lengths of type elements
  MPI_Datatype block_types [] =
  {
    MPI_UNSIGNED_CHAR,
    MPI_UNSIGNED_CHAR,
    MPI_UNSIGNED_CHAR
  }; //Set types
  MPI_Aint start, displ[3];

  MPI_Address( &pixel, &start );
  MPI_Address( &pixel->r, &displ[0] );
  MPI_Address( &pixel->g, &displ[1] );
  MPI_Address( &pixel->b, &displ[2] );

  displ[0] -= start; // Displacement relative to address of start
  displ[1] -= start; // Displacement relative to address of start
  displ[2] -= start; // Displacement relative to address of start

  MPI_Type_struct( 3, block_lengths, displ, block_types, &mpi_pixel );
  MPI_Type_commit( &mpi_pixel );
}
