#include "read_file.hpp"

using namespace std;

size_data read_ppm_width_height(string filename)
{
  FILE* file_stream;
  file_stream = fopen(filename.c_str(),"r");

  if (file_stream==NULL) {
    perror ("Error opening file");
    exit(1);
  }

  size_data data;
  data.width = ppm_readint(file_stream);
  data.height = ppm_readint(file_stream);

  return data;
}

// -------------
