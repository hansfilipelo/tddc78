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
        // Check for collisions.
        // Move particles that has not collided with another.
        // Check for wall interaction and add the momentum.
        // Check if particles need to be communicated
        pcord_t* particle;
        pcord_t* other_particle;
        size_t last_pos = particles.size()-1;
        bool has_collision = false;

        if (particles.size() > 0){
            for (size_t i = 0; i < last_pos; i++) {
                has_collision = false;

                particle = particles.at(i);
                for (size_t j = i+1; j < last_pos+1; j++) {

                    other_particle = particles.at(j);
                    collision = collide(particle, other_particle);

                    if (collision == -1) {
                        continue;
                    }
                    else{
                        interact(particle, other_particle, collision);

                        total_momentum += wall_collide(particle,box);
                        total_momentum += wall_collide(other_particle,box);

                        tmp_particles.push_back(other_particle);
                        tmp_particles.push_back(particle);

                        swap(other_particle, particles.at(last_pos));
                        last_pos--;
                        has_collision = true;
                        break;
                    }
                } // End loop with iterator j

                if (has_collision == false) {
                    if ( particle->y > my_cords.y0 && particle->y < my_cords.y1 ) {
                        feuler(particle, STEP_SIZE);
                        tmp_particles.push_back(particle);
                        total_momentum += wall_collide(particle,box);
                    }
                    else if ( particle->y < my_cords.y0 ){
                        up_transfers.push_back(*particle);
                        pos_to_erase.push_back(i);
                    }
                    else{
                        down_transfers.push_back(*particle);
                        pos_to_erase.push_back(i);
                    }
                }
            } // End loop with iterator i

            // Check last particle as well. It can not collide with anything.
            particle = particles.at(last_pos);
            if (has_collision == false) {
                if ( particle->y > my_cords.y0 && particle->y < my_cords.y1 ) {
                    feuler(particle, STEP_SIZE);
                    tmp_particles.push_back(particle);
                    total_momentum += wall_collide(particle,box);
                }
                else if ( particle->y < my_cords.y0 ){
                    up_transfers.push_back(*particle);
                    pos_to_erase.push_back(last_pos);
                }
                else{
                    down_transfers.push_back(*particle);
                    pos_to_erase.push_back(last_pos);
                }
            }

            // Free memory
            size_t pos_to_erase_size = pos_to_erase.size();

            for (size_t i = 0; i < pos_to_erase_size; i++) {
                delete particles.at(pos_to_erase.back());
                pos_to_erase.pop_back();
            }
            particles.erase(particles.begin(), particles.end());
            particles.swap(tmp_particles);
        }

        // Send particles who should change processing element
        if (my_rank != n_tasks-1) {
            send_count = down_transfers.size();
            MPI_Isend(&send_count, 1, MPI_UNSIGNED, my_rank+1, 2*(my_rank+1), com, &send_count_request);

            if(send_count != 0) {
                MPI_Isend(&down_transfers.at(0), send_count, mpi_particle, my_rank+1, my_rank+1, com, &send_data_request);
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
                size_t transfer_size;
                vector<int> back_to_particles;
                last_pos = up_transfers.size() - 1;

                for (size_t i = 0; i < recv_count; i++) {
                    particle = &recv_buffer[i];

                    for (size_t j = 0; j < last_pos+1; j++) {
                        other_particle = &up_transfers.at(j);

                        collision = collide(particle, other_particle);

                        if (collision != -1) {
                            interact(particle, other_particle, collision);

                            total_momentum += wall_collide(particle,box);
                            total_momentum += wall_collide(other_particle,box);

                            if ( other_particle->y > my_cords.y0 && other_particle->y < my_cords.y1 ) {
                                particles.push_back(Utils::copy_particle(*other_particle));
                                Utils::pcord_swap(other_particle, &up_transfers.at(last_pos));
                                up_transfers.erase(up_transfers.begin()+last_pos);
                            }
                            else {
                                Utils::pcord_swap(other_particle, &up_transfers.at(last_pos));
                            }
                            last_pos--;

                            if ( particle->y > my_cords.y0 && particle->y < my_cords.y1 ){
                                particles.push_back(Utils::copy_particle(*particle));
                            }
                            else{
                                up_transfers.push_back(*particle);
                            }

                            has_collision = true;
                            break;
                        }
                    } // End of loop with iterator j

                    if (has_collision == false) {
                        feuler(particle, STEP_SIZE);
                        total_momentum += wall_collide(particle,box);
                        particles.push_back(Utils::copy_particle(*particle));
                    }
                } // End of loop with iterator i

                for (size_t i = 0; i < last_pos+1; i++) {
                    particle = &up_transfers.at(i);
                    feuler(particle, STEP_SIZE);
                    total_momentum += wall_collide(particle, box);
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
                MPI_Isend(&up_transfers.at(0), send_count, mpi_particle, my_rank-1, my_rank-1, com, &send_data_request);
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
                    particles.push_back(Utils::copy_particle(recv_buffer[i]));
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
