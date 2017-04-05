/*SUBROUTINE SAVE_TO_FILE(array,xdim,ydim,zdim,name,name_length)

  **************************************************************

  This subroutine saves images as binary files with NO HEADER.
  Each voxel value in output image is a
  4 BYTE FLOATING POINT number.

  **************************************************************
  This is simply a subroutine which saves the phantom to a file.
  This subroutine can be modified by the user to save the
  phantom any whatever format the user wants. The purpose of
  this seperate subroutine is to allow various users to
  modify the file format of the saved phantom without having
  to modify any of the other files containing the phantom code.
  **************************************************************
	
  array:           contains the data to be saved
  xdim,ydim,zdim:  dimensions of data
  name:            name of file
  --------------------------------------------------------------
*/
/* STANDARD C LIBRARIES */
#include "global_includes.h"
#include <string>

typedef float FLOATTYPE;

/* Not currently used, but here if you need them: */
/* typedef unsigned char BYTETYPE;   /* value range: 0 to 255 */
/* typedef short INTTYPE;            /* value range: -32768 to 32767 */

/*-------------------------------------------------------------*/
void SAVE_TO_FILE(float *array, int xdim, int ydim, int zdim, char *name)
{
  int	n, TotalPix;
  char  new_name[64];
  FILE	*fp_out;
/*
** ADD EXTENSION TO FILENAME:  
*/
  strcpy(new_name, name);
  strcat(new_name, ".bin");   /* add extension to file name*/
  unlink(new_name);          
/*
** WRITE TO FILE: 
*/
  TotalPix = xdim*ydim*zdim;   /* total number of pixels in 3D image */
  if ((fp_out = fopen(new_name, "wb")) == NULL)
    Abort("Cannot open raw output file");
   
  n = fwrite (array, sizeof(FLOATTYPE), TotalPix, fp_out);
  if (n != TotalPix) {
    printf("Error : fwrite return %d\n",n);
    Abort ("Failure writing pixels to output image\n");
    unlink (new_name);
  }
  fclose(fp_out);
}
