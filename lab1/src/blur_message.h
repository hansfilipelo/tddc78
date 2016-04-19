#include <mpi.h>


typedef struct {
  unsigned int radius;
  unsigned int width;
  unsigned int height;
  char* data[MAX_PIXELS*3];
} blur_message;
