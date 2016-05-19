#include "Utils.hpp"

using namespace std;

void Utils::init(int my_rank)
{
    srand(time(0)+my_rank);
}

float Utils::generate_random_float(float from, float to)
{
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = to - from;
    float r = random * diff;
    return from + r;
}

pcord_t Utils::init_particle(cord_t my_cords)
{
    pcord_t p;
    p.x = Utils::generate_random_float(my_cords.x0, my_cords.x1);
    p.y = Utils::generate_random_float(my_cords.y0, my_cords.y1);
    p.vx = Utils::generate_random_float(- MAX_INITIAL_VELOCITY, MAX_INITIAL_VELOCITY);
    p.vy = Utils::generate_random_float(- MAX_INITIAL_VELOCITY, MAX_INITIAL_VELOCITY);

    return p;
}

void Utils::pcord_swap(pcord_t* p1, pcord_t* p2)
{
    pcord_t tmp;

    tmp.x = p1->x;
    tmp.y = p1->y;
    tmp.vx = p1->vx;
    tmp.vy = p1->vy;

    p1->x = p2->x;
    p1->y = p2->y;
    p1->vx = p2->vx;
    p1->vy = p2->vy;

    p2->x = tmp.x;
    p2->y = tmp.y;
    p2->vx = tmp.vx;
    p2->vy = tmp.vy;
}
