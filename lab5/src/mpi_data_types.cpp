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

  displ[0] = 0;
  displ[1] = sizeof(float);
  displ[2] = 2*sizeof(float);
  displ[3] = 3*sizeof(float);

  MPI_Type_create_struct( 4, block_lengths, displ, block_types, mpi_particle );
  MPI_Type_commit( mpi_particle );
}
