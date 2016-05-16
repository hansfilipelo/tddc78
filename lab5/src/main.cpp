#include "main.hpp"

using namespace std;

int main(int argc, char** argv)
{
    // Init mpi
    int n_tasks, my_rank;
    MPI_Datatype mpi_particle;
    pcord_t* mpi_p;

    MPI_Request send_count_request;
    MPI_Request send_data_request;
    MPI_Request receive_data_request;
    MPI_Request receive_count_request;
    unsigned recv_count = 0;
    int send_count = 0;
    pcord_t* recv_buffer;

    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(com, &n_tasks);
    MPI_Comm_rank(com, &my_rank);
    create_mpi_particle_t(mpi_p, &mpi_particle);

    // Init random nr generator
    Utils::init(my_rank);

    // Initiate walls for box as well as split box into sub-areas
    const cord_t box = {0, BOX_HORIZ_SIZE, 0, BOX_VERT_SIZE};
    vector<pcord_t*> particles;
    vector<pcord_t*> tmp_particles;
    float vert_stop;

    if ( my_rank != n_tasks-1 ) {
        vert_stop = ((float)BOX_VERT_SIZE/n_tasks)*(my_rank+1);
    }
    else {
        vert_stop = (float)BOX_VERT_SIZE;
    }

    cord_t my_cords = {0, (float)BOX_HORIZ_SIZE, ((float)BOX_VERT_SIZE/n_tasks)*my_rank, vert_stop};

    // Initiate particles
    for (size_t i = 0; i < INIT_NO_PARTICLES; i++) {
        pcord_t* p = new pcord_t();
        p->x = Utils::generate_random_float(my_cords.x0, my_cords.x1);
        p->y = Utils::generate_random_float(my_cords.y0, my_cords.y1);
        p->vx = Utils::generate_random_float(0, MAX_INITIAL_VELOCITY);
        p->vy = Utils::generate_random_float(0, MAX_INITIAL_VELOCITY);
        particles.push_back(p);
    }

    // Main loop: for each time-step do
    float total_momentum = 0;
    float collision;

    // Temporary storage for transfers
    vector<pcord_t> up_transfers;
    vector<pcord_t> down_transfers;

    vector<int> pos_to_erase;

    for (size_t t = 0; t < _SIMULATION_STEPS_; t++) {
        cout << "My rank: " << my_rank << ", sim step: " << t << endl;
        // Check for collisions.
        // Move particles that has not collided with another.
        // Check for wall interaction and add the momentum.
        // Check if particles need to be communicated
        pcord_t* particle;
        pcord_t* other_particle;
        size_t last_pos = particles.size()-1;

        for (size_t i = 0; i < last_pos; i++) {
            particle = particles.at(i);
            for (size_t j = i+1; j < last_pos+1; j++) {
                other_particle = particles.at(j);
                collision = collide(particle, other_particle);
                interact(particle, other_particle, collision);
            }

            total_momentum += wall_collide(particle,box);

            if ( particle->y < my_cords.y0 ) {
                up_transfers.push_back(*particle);
                pos_to_erase.push_back(i);
            }
            else if ( particle->y > my_cords.y1 ){
                down_transfers.push_back(*particle);
                pos_to_erase.push_back(i);
            }
            else{
                tmp_particles.push_back(particle);
            }
        }

        particle = particles.at(last_pos);
        total_momentum += wall_collide(particle,box);

        if ( particle->y < my_cords.y0 ) {
            up_transfers.push_back(*particle);
            pos_to_erase.push_back(last_pos);
        }
        else if ( particle->y > my_cords.y1 ){
            down_transfers.push_back(*particle);
            pos_to_erase.push_back(last_pos);
        }
        else{
            tmp_particles.push_back(particle);
        }

        // Free memory
        size_t pos_to_erase_size = pos_to_erase.size();
        for (size_t i = 0; i < pos_to_erase_size; i++) {
            delete particles.at(pos_to_erase.back());
            pos_to_erase.pop_back();
        }
        particles.erase(particles.begin(), particles.end());
        particles.swap(tmp_particles);

        // Send particles who should change processing element
        if (my_rank != n_tasks-1) {
            send_count = down_transfers.size();
            MPI_Isend(&send_count, 1, MPI_UNSIGNED, my_rank+1, 2*(my_rank+1), com, &send_count_request);

            if(send_count != 0) {
                cout << "Before first ibsend for rank " << my_rank << endl;
                MPI_Isend(&down_transfers.at(0), send_count, mpi_particle, my_rank+1, my_rank+1, com, &send_data_request);
                cout << "After first ibsend for rank " << my_rank << endl;
            }
        }

        if ( my_rank != 0) {
            // Receive the nr of elements to get from processor my_rank-1
            MPI_Irecv(&recv_count, 1, MPI_UNSIGNED, my_rank-1, 2*my_rank, com, &receive_count_request);
            MPI_Wait(&receive_count_request, MPI_STATUS_IGNORE);

            if(recv_count != 0) {
                // Allocate buffer
                recv_buffer = (pcord_t*)malloc(sizeof(pcord_t)*recv_count);

                // Receive elements from my_rank-1
                MPI_Irecv(recv_buffer, recv_count, mpi_particle, my_rank-1, my_rank, com, &receive_data_request);
                MPI_Wait(&receive_data_request, MPI_STATUS_IGNORE);

                // Check whether receive particles and upgoing particles collide
                pcord_t *particle, *other_particle;
                size_t transfer_size = up_transfers.size();
                vector<int> back_to_particles;

                for (size_t i = 0; i < recv_count; i++) {
                    for (size_t j = 0; j < transfer_size; j++) {

                        particle = &recv_buffer[i];
                        other_particle = &up_transfers.at(j);

                        collision = collide(particle, other_particle);
                        interact(particle, other_particle, collision);

                        // Check if incomming particle should be sent back or kept
                        if ( particle->y < my_cords.y0 ) {
                            up_transfers.push_back(*particle);
                        }
                        else{
                            pcord_t* p = new pcord_t();
                            p->x = particle->x;
                            p->y = particle->y;
                            p->vx = particle->vx;
                            p->vy = particle->vy;
                            particles.push_back(p);
                        }
                        // Check if we should keep outgoing particle
                        if (other_particle->y >= my_cords.y0) {
                            back_to_particles.push_back(j);
                        }
                    }
                }

                // Put back outgoing particles if they should be put back in particles
                transfer_size = back_to_particles.size();
                unsigned pos;

                for (size_t i = 0; i < transfer_size; i++) {
                    pos = back_to_particles.back();
                    pcord_t* p = new pcord_t();
                    p->x = recv_buffer[pos].x;
                    p->y = recv_buffer[pos].y;
                    p->vx = recv_buffer[pos].vx;
                    p->vy = recv_buffer[pos].vy;
                    particles.push_back(p);
                    back_to_particles.pop_back();
                }

                // Free the receive buffer
                free(recv_buffer);
            }


            // Wait for send to finish if not finished
            if (my_rank != n_tasks-1) { // First see that send_count is available for reuse
                MPI_Wait(&send_count_request, MPI_STATUS_IGNORE);
                if(send_count != 0) {
                    MPI_Wait(&send_data_request, MPI_STATUS_IGNORE);
                }
            }

            send_count = up_transfers.size();
            MPI_Isend(&send_count, 1, MPI_UNSIGNED, my_rank-1, 2*my_rank, com, &send_count_request);
            if(send_count != 0) {
                cout << "Before second ibsend for rank " << my_rank << ", send count: " << send_count << endl;
                MPI_Isend(&up_transfers.at(0), send_count, mpi_particle, my_rank-1, my_rank-1, com, &send_data_request);
                cout << "After second ibsend for rank " << my_rank << endl;
            }
        }


        // Receive particles from my_rank+1
        if (my_rank != n_tasks-1) {
            // Receive the nr of elements to get from processor my_rank+1
            MPI_Irecv(&recv_count, 1, MPI_UNSIGNED, my_rank+1, 2*(my_rank+1), com, &receive_count_request);
            MPI_Wait(&receive_count_request, MPI_STATUS_IGNORE);

            if(recv_count != 0) {
                // Allocate recv buffer
                recv_buffer = (pcord_t*)malloc(sizeof(pcord_t)*recv_count);

                // Get the particles!
                MPI_Irecv(recv_buffer, recv_count, mpi_particle, my_rank+1, my_rank, com, &receive_data_request);
                MPI_Wait(&receive_data_request, MPI_STATUS_IGNORE);

                for (size_t i = 0; i < recv_count; i++) {
                    pcord_t* p = new pcord_t();
                    p->x = recv_buffer[i].x;
                    p->y = recv_buffer[i].y;
                    p->vx = recv_buffer[i].vx;
                    p->vy = recv_buffer[i].vy;
                    particles.push_back(p);
                }

                free(recv_buffer);
            }
        }

        // Wait for send stuff to get free for next run.
        MPI_Wait(&send_count_request, MPI_STATUS_IGNORE);
        if(send_count != 0) {
            MPI_Wait(&send_data_request, MPI_STATUS_IGNORE);
        }
        up_transfers.clear();
        down_transfers.clear();
    }

    // Reduction and calculate pressure.
    total_momentum = total_momentum/(_SIMULATION_STEPS_ * STEP_SIZE);

    if(my_rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, &total_momentum, 1, MPI_FLOAT, MPI_SUM, 0, com);
            cout << "Pressure = " << total_momentum << endl;
    }
    else {
        MPI_Reduce(&total_momentum, &total_momentum, 1, MPI_FLOAT, MPI_SUM, 0, com);
    }

    MPI_Finalize();
    return 0;
}
