#include "Utils.hpp"

using namespace std;

void Utils::init()
{
    srand(time(0));
}

float Utils::generate_random_float(float from, float to)
{
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = to - from;
    float r = random * diff;
    return from + r;
}
