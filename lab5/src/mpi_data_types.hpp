#ifndef _MPI_DATA_TYPES_H_
#define _MPI_DATA_TYPES_H_

#include <mpi.h>
#include "definitions.hpp"

void create_mpi_particle_t(const pcord_t* particle, MPI_Datatype* mpi_particle);

#endif // _MPI_DATA_TYPES_H_
