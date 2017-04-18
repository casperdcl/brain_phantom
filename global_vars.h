/*---------------------------------------------
**       GLOBAL VARIABLES
**---------------------------------------------*/

#pragma once

#define q_degree 3
#define p_degree 3

extern const float PI;
extern int gender_flag;

/* GENERAL PARAMETERS */
extern int 	act_phan_ave;
extern int 	act_phan_each;
extern int 	motion_option;
extern float 	out_period;
extern int 	out_frames;
extern float 	motion_period;
extern float	motion_start_ph_index;
extern float	x_trans, y_trans, z_trans;
extern float	x_rot, y_rot, z_rot;

extern float  	pixel_width;
extern float  	slice_width;
extern int 	array_size;
extern int 	subvoxel_index;
extern int 	subvxl_index;
extern int 	startslice;
extern int 	endslice;

/*Array of surfaces for torso models*/
extern TRI_MODEL tmodel[256];

/* For pixelizing the NURBS surfaces */
extern int   txdim, tydim, tzdim;
extern int   txdim2, tydim2, tzdim2;
extern int   xoffset, yoffset, zoffset; /*Used to center individual structures for temporary rendering*/
extern int   xoffset2, yoffset2, zoffset2; /*Used to center individual structures for temporary rendering*/
extern float xspan, yspan, zspan;

extern float   xoff, yoff, zoff; /*Used to center torso in the final 3D image produced*/
extern int     XDIM,YDIM,ZDIM,ZDIM_OUTPUT,slice_index;

extern float voxel_volume,subvxl_vol;

extern float activity[256];
extern int motion_flag;
extern int user_time_flag;
extern char user_motion_file[200];
extern CURVE motion_curve[6];
