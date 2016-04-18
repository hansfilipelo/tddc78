#include "ppm.hpp"

using namespace std;

void ppm_error(const char *str)
{
  fprintf(stderr, "%s\n", str);
  exit(1);
}

/* Read a character from a ppm-file. Remove all comments */
char ppm_readchar(FILE *file)
{
  string error_message = "read error";

  char ch;
  ch = getc(file);
  if (ch==EOF)
  ppm_error(error_message.c_str());
  if (ch=='#') do {
    ch = getc(file);
    if (ch==EOF)
    ppm_error(error_message.c_str());
  } while (ch != '\n');

  return ch;
}

/* Read the magic number in a ppm-file */
int ppm_readmagicnumber(FILE *file)
{
  int ch1, ch2;
  string error_message = "read error";

  ch1 = getc(file);
  if (ch1==EOF)
  ppm_error(error_message.c_str());
  ch2 = getc(file);
  if (ch2==EOF)
  ppm_error(error_message.c_str());
  return ch1 * 256 + ch2;
}

/* Read an ASCII integer from a ppm-file */
int ppm_readint(FILE *file)
{
  char ch;
  int i;
  string error_message = "error in readint";

  do
  ch = ppm_readchar(file);
  while (ch == ' ' || ch == '\t' || ch == '\n');

  if (ch < '0' || ch > '9')
  ppm_error(error_message.c_str());
  i = 0;
  do {
    i = i*10 + (ch - '0');
    ch = ppm_readchar(file);
  } while (ch >= '0' && ch <= '9');
  return i;
}
