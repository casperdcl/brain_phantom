/*	SUBROUTINE SAVE_TO_FILE(array,xdim,ydim,zdim,name,name_length)

	**************************************************************

	This file saves image in UNC's .IM format
	The extension .im is added to the file name

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
#include "image.h"

#define Abort(Mesg){fprintf(stderr, "%s\n", Mesg); exit(1);}

void	SAVE_TO_FILE(float *array, int xdim, int ydim, int zdim, char *name)
{
	char		new_name[64];
	int 		dimv[3];
	int		dimc = 3;
	IMAGE		*out_image;

	dimv[0] = zdim;
	dimv[1] = ydim;
	dimv[2] = xdim;

	
/*	ADD EXTENSION TO FILENAME:	*/
	strcpy(new_name, name);
	strcat(new_name, ".im");   /*extension*/
        unlink(new_name);

/*	SAVE FILE:	*/

        /* Create the output image file */
        if ((out_image = imcreat(new_name, DEFAULT, REAL, dimc, dimv)) == INVALID)
             Abort("Can not create output image");

        /* Write out pixel values to the new image file */
        if (imwrite(out_image, 0, xdim*ydim*zdim-1, array) == INVALID)
           Abort("Can not write pixels to output image\n");

	imclose(out_image);

}
