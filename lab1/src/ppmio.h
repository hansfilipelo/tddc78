/*
  File: ppmio.h

  Declarations for the PPM image file IO functions.

*/
#ifndef _PPMIO_H_
#define _PPMIO_H_

/* Function: read_ppm - reads data from an image file in PPM format.
   Input: fname - name of an image file in PPM format to read.
   Output:
      xpix, ypix - size of the image in x & y directions
      max - maximum intensity in the picture
      data - color data array. MUST BE PREALLOCATED to at least MAX_PIXELS*3 bytes.
   Returns: 0 on success.
 */
int read_ppm (const char * fname,
	       unsigned int * xpix, unsigned int * ypix, int * max, char * data);

/* Function: write_ppm - write out an image file in PPM format.
   Input:
      fname - name of an image file in PPM format to write.
      xpix, ypix - size of the image in x & y directions
      data - color data.
   Returns: 0 on success.
 */
int write_ppm (const char * fname, unsigned int xpix, unsigned int ypix, char * data);

#endif
