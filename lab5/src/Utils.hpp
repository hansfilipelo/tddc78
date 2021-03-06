#include <ctime>
#include <cstdlib>
#include <iostream>
#include "coordinate.hpp"
#include "definitions.hpp"

#define UP true
#define DOWN false

class Utils
{
public:
    static float generate_random_float(float from, float to);
    static void init(int my_rank);
    static pcord_t init_particle(cord_t my_cords);
    static void pcord_swap(pcord_t* p1, pcord_t* p2);
    static bool will_pass_edge(pcord_t* particle, float edge, bool transfer);
private:
};
