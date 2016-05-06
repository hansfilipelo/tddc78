#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"

#ifdef __APPLE__
#include "timing_mach.h"
#include "barrier.h"
#endif

#ifdef _linux_
#include <sys/time.h>
#endif

// ----------

int radius, colmax;
int offset;
double* w;

struct timespec stime, etime;

size_data_t size_data;
unsigned int partitioned_height;
unsigned int remainder_height;
pixel_t* src;


// Memory barrier
pthread_barrier_t x_done_barrier;

pthread_t worker_handles[_N_TASKS_];

// ----------

void* pthread_blur_filter(void* id)
{
    int my_id = *(int*)id;

    blurfilter(size_data.width, size_data.height, partitioned_height, src, radius, w, my_id, _N_TASKS_, &x_done_barrier);

    exit(0);
}


// ----------

int main (int argc, char ** argv) {
    // Init memory barrier
    int i; // iterator for later
    pthread_barrier_init(&x_done_barrier,NULL,_N_TASKS_);

    /* Take care of the arguments */
    if (argc != 4) {
        fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
        exit(1);
    }

    radius = atoi(argv[1]);

    if((radius > MAX_RAD) || (radius < 1)) {
        fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
        exit(1);
    }

    w = malloc(sizeof(double)*radius);

    /* Read input image */
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

    // Calculate some sweet height values
    partitioned_height = size_data.height/_N_TASKS_;
    remainder_height = size_data.height-(_N_TASKS_*partitioned_height);

    /* Calculate gaussian weights */
    get_gauss_weights(radius, w);

    printf("Start worker threads to filter image\n");
    clock_gettime(CLOCK_REALTIME, &stime);


    // Spawn worker threads running the filter
    int thread_ids[_N_TASKS_];
    for (i = 0; i < _N_TASKS_; i++) {
        thread_ids[i] = i;
        pthread_create(&worker_handles[i], NULL, pthread_blur_filter, &thread_ids[i]);
    }

    // Wait for threads to finnish
    for (i = 0; i < _N_TASKS_; i++) {
        pthread_join(worker_handles[i], NULL);
    }

    // Time
    clock_gettime(CLOCK_REALTIME, &etime);

    printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
    1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");

    if(write_ppm (argv[3], size_data.width, size_data.height, (char *)src) != 0){
        exit(1);
    }

    return(0);
}
