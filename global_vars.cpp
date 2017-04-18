/*---------------------------------------------
**       GLOBAL VARIABLES
**---------------------------------------------*/

#include "nurbs.h"

extern const float PI = 3.14159265358979323846264338328f;
int gender_flag;

/* GENERAL PARAMETERS */
int 	act_phan_ave;
int 	act_phan_each;
int 	motion_option;
float 	out_period;
int 	out_frames;
float 	motion_period;
float	motion_start_ph_index;
float	x_trans, y_trans, z_trans;
float	x_rot, y_rot, z_rot;

float  	pixel_width;
float  	slice_width;
int 	array_size;
int 	subvoxel_index;
int 	subvxl_index;
int 	startslice;
int 	endslice;

/*Array of surfaces for torso models*/
TRI_MODEL tmodel[256];

/* For pixelizing the NURBS surfaces */
int   txdim, tydim, tzdim;
int   txdim2, tydim2, tzdim2;
int   xoffset, yoffset, zoffset; /*Used to center individual structures for temporary rendering*/
int   xoffset2, yoffset2, zoffset2; /*Used to center individual structures for temporary rendering*/
float xspan, yspan, zspan;

float   xoff, yoff, zoff; /*Used to center torso in the final 3D image produced*/
int     XDIM,YDIM,ZDIM,ZDIM_OUTPUT,slice_index;

float voxel_volume,subvxl_vol;

float activity[256];
int motion_flag = 0;
int user_time_flag = 0;
char user_motion_file[200];
CURVE motion_curve[6];
