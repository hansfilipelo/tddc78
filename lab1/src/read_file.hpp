#include <iostream>
#include <string>
#include <vector>
#include <cstdio>

#include "ppm.hpp"
#include "ppmio.hpp"

struct size_data {
  unsigned int width;
  unsigned int height;
};

size_data read_ppm_width_height(std::string filename);


char* read_ppm_file(std::string filename, unsigned int width, unsigned int height);
