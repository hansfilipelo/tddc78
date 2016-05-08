#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#include "ppmio.h"
#include "thresfilter.h"
#include "gaussw.h"

#ifdef __APPLE__
#include "timing_mach.h"
#include "barrier.h"
#endif

#ifdef _linux_
#include <sys/time.h>
#endif

// ----------

int colmax;
int offset;

struct timespec start_time, end_time;

size_data_t size_data;
unsigned int total_pixels;
unsigned int partitioned_height;
unsigned int remainder_height;
pixel_t* src;
double intermediate_means[_N_TASKS_+1];
double total_mean = -1;


// Memory barrier
pthread_barrier_t int_mean_done_barrier;
pthread_barrier_t total_mean_done_barrier;
pthread_barrier_t thres_done_barrier;

pthread_t worker_handles[_N_TASKS_];

// ----------

void* pthread_thres_filter(void* id)
{
    int my_id = *(int*)id;

    intermediate_means[my_id] = calculate_partitioned_mean(partitioned_height*size_data.width, total_pixels, src, my_id);

    pthread_barrier_wait(&int_mean_done_barrier);

    pthread_barrier_wait(&total_mean_done_barrier);

    thresfilter(partitioned_height*size_data.width, total_pixels, total_mean, src, my_id);

    pthread_barrier_wait(&thres_done_barrier);

    return NULL;
}


// ----------

int main (int argc, char ** argv) {
    // Init memory barrier
    int i; // iterator for later
    pthread_barrier_init(&int_mean_done_barrier,NULL,_N_TASKS_+1);
    pthread_barrier_init(&total_mean_done_barrier,NULL,_N_TASKS_+1);
    pthread_barrier_init(&thres_done_barrier, NULL, _N_TASKS_+1);

    /* Take care of the arguments */
    if (argc != 3) {
        fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
        exit(1);
    }

    /* Read input image */
    src = malloc(sizeof(pixel_t) * MAX_PIXELS);

    /* read file */
    if(read_ppm (argv[1], &size_data.width, &size_data.height, &colmax, (char *) src) != 0) {
        exit(1); // TODO: Broadcast exit
    }
    total_pixels = size_data.width*size_data.height;

    if (colmax > 255) {
        fprintf(stderr, "Too large maximum color-component value\n");
        exit(1); // TODO: Broadcast exit
    }

    printf("Has read the image, generating coefficients\n");

    // Calculate some sweet height values
    partitioned_height = size_data.height/_N_TASKS_;
    remainder_height = size_data.height-(_N_TASKS_*partitioned_height);

    printf("Start worker threads to filter image\n");
    clock_gettime(CLOCK_REALTIME, &start_time);


    // Spawn worker threads running the filter
    int thread_ids[_N_TASKS_];
    for (i = 0; i < _N_TASKS_; i++) {
        thread_ids[i] = i;
        pthread_create(&worker_handles[i], NULL, pthread_thres_filter, &thread_ids[i]);
    }

    intermediate_means[_N_TASKS_] = calculate_partitioned_mean(partitioned_height*size_data.width, total_pixels, src, _N_TASKS_);

    // Wait for threads to finnish
    pthread_barrier_wait(&int_mean_done_barrier);

    printf("Here 1.\n");

    total_mean = calculate_total_mean(intermediate_means, _N_TASKS_+1, partitioned_height*size_data.width, remainder_height*size_data.width, total_pixels);

    pthread_barrier_wait(&total_mean_done_barrier);

    printf("Here 2.\n");

    pthread_barrier_wait(&thres_done_barrier);

    printf("Here 3.\n");

    // Time
    clock_gettime(CLOCK_REALTIME, &end_time);

    printf("Filtering took: %g secs\n", (end_time.tv_sec  - start_time.tv_sec) +
    1e-9*(end_time.tv_nsec  - start_time.tv_nsec)) ;

    /* write result */
    printf("Writing output file\n");

    if(write_ppm (argv[2], size_data.width, size_data.height, (char *)src) != 0){
        exit(1);
    }

    return(0);
}
