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

pcord_t* Utils::copy_particle(pcord_t particle)
{
    pcord_t* p = new pcord_t();
    p->x = particle.x;
    p->y = particle.y;
    p->vx = particle.vx;
    p->vy = particle.vy;

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
