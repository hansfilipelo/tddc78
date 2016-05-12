#include "main.hpp"

using namespace std;

int main(int argc, char** argv)
{
    // Init random nr generator
    Utils::init();
    // Init mpi
    int n_tasks, my_rank;
    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(com, &n_tasks);
    MPI_Comm_rank(com, &my_rank);

    // Initiate walls for box as well as split box into sub-areas
    const cord_t box = {0, BOX_HORIZ_SIZE, 0, BOX_VERT_SIZE};
    vector<particle_t> particles;
    float vert_stop;

    if ( my_rank != n_tasks-1 ) {
        vert_stop = ((float)BOX_VERT_SIZE/n_tasks)*(my_rank+1);
    }
    else {
        vert_stop = (float)BOX_VERT_SIZE;
    }

    cord_t my_cords = {0, (float)BOX_HORIZ_SIZE, ((float)BOX_HORIZ_SIZE/n_tasks)*my_rank, vert_stop};

    // Initiate particles
    for (size_t i = 0; i < INIT_NO_PARTICLES; i++) {
        particle_t* p = new particle_t();
        p->pcord.x = Utils::generate_random_float(my_cords.x0, my_cords.x1);
        p->pcord.y = Utils::generate_random_float(my_cords.y0, my_cords.y1);
        p->pcord.vx = Utils::generate_random_float(0, MAX_INITIAL_VELOCITY);
        p->pcord.vy = Utils::generate_random_float(0, MAX_INITIAL_VELOCITY);
        particles.push_back(*p);
    }

    // Main loop: for each time-step do
    float total_momentum = 0;
    float collision;
    for (size_t t = 0; t < _SIMULATION_STEPS_; t++) {

        for (vector<particle_t>::iterator particle = particles.begin(); particle != particles.end()-1; ++particle) {

            for (vector<particle_t>::iterator other_particle = particle+1; other_particle != particles.end(); ++other_particle) {

                collision = collide(&particle->pcord, &other_particle->pcord);
                interact(&particle->pcord, &other_particle->pcord, collision);
            }

            total_momentum += wall_collide(&particle->pcord,box);
        }
    }

    printf("%f\n", total_momentum);
    // for all particles do
    // Check for collisions.
    // Move particles that has not collided with another.
    // Check for wall interaction and add the momentum.
    // Communicate if needed.
    // Calculate pressure.

    MPI_Finalize();
    return 0;
}
