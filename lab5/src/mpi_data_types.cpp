#include "mpi_data_types.hpp"

// Create mpi size data type
void create_mpi_particle_t(const pcord_t* particle, MPI_Datatype* mpi_particle)
{
  int block_lengths [] = { 1, 1, 1, 1 }; // Lengths of type elements
  MPI_Datatype block_types [] =
  {
    MPI_FLOAT,
    MPI_FLOAT,
    MPI_FLOAT,
    MPI_FLOAT
  }; //Set types
  MPI_Aint start, displ[4];

  MPI_Get_address( particle, &start );
  MPI_Get_address( &particle->x, &displ[0] );
  MPI_Get_address( &particle->y, &displ[1] );
  MPI_Get_address( &particle->vx, &displ[2] );
  MPI_Get_address( &particle->vy, &displ[3] );

  displ[0] -= start; // Displacement relative to address of start
  displ[1] -= start; // Displacement relative to address of start
  displ[2] -= start;
  displ[3] -= start;

  MPI_Type_create_struct( 4, block_lengths, displ, block_types, mpi_particle );
  MPI_Type_commit( mpi_particle );
}
