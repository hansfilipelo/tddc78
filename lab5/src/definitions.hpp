#include<stdlib.h>
#include<math.h>

#include "coordinate.hpp"
#include "physics.hpp"

#ifndef _definitions_h
#define _definitions_h

#define PI 3.141592653

#define MAX_NO_PARTICLES  40000  /* Maximum number of particles/processor */
#define INIT_NO_PARTICLES 10000   /* Initial number of particles/processor */
#define MAX_INITIAL_VELOCITY 50

#define BOX_HORIZ_SIZE 10000.0
#define BOX_VERT_SIZE 10000.0
#define WALL_LENGTH (2.0*BOX_HORIZ_SIZE+2.0*BOX_VERT_SIZE)

#define PARTICLE_BUFFER_SIZE MAX_NO_PARTICLES/5
#define COMM_BUFFER_SIZE  5*PARTICLE_BUFFER_SIZE

typedef struct {
  pcord_t  pcord;
  int ptype;        /* Used to simulate mixing of gases */
} particle_t;



#endif
