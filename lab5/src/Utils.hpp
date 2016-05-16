#include <ctime>
#include <cstdlib>
#include "coordinate.hpp"

class Utils
{
public:
    static float generate_random_float(float from, float to);
    static void init(int my_rank);
    static pcord_t* copy_particle(pcord_t particle);
    static void pcord_swap(pcord_t* p1, pcord_t* p2);
private:
};
