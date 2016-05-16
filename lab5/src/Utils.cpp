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
