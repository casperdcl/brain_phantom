static const char *const DOCUMENTATION =
R"(Program to produce voxelised phantoms from NURBS torso file input

Usage:
  %s [options] <gen_par> <out_base>

Options:
  -x X  translation in mm
  -y Y  translation in mm
  -z Z  translation in mm
  -X X  rotation in degrees
  -Y Y  rotation in degrees
  -Z Z  rotation in degrees

Arguments:
  <gen_par>   : general parameter filename. Note that start_slice should be
                an even number > 1 to avoid artefacts.
  <out_base>  : base string for output files

Paul Segars <paul.segars@duke.edu> 2006-09-03
Casper da Costa-Luis <casper.dcl@kcl.ac.uk> 2017-04
)";

#include "global_includes.h"
#include "nurbs.h"            /* header file for nurbs */
#include "global_vars.h"      /* global variables for ncat */
#include "constants.h"        /* global constants */
#include "dncatsubs.h"
#include "dncat_output.h"

#include <cstring>  // strcpy, strcat, memcpy

char  gen_parfile[64];
char  outputbase[64];

float    *fphan, *fphan_id, *fphan_avg, *slice;
TRI_MODEL slice_tmodel;

/***********************/
/* BEGIN MAIN PROGRAM  */
/***********************/

int main(int argc, char *argv[])
{
  setbuf(stdout, NULL);

  FILE *fp, *fp2;
  int i, j, k, l, kf;
  char filename[200], outfile[200], logfile[200];
  unsigned long index1, index2;
  int intensity;
  float centerx, centery, centerz;
  float XCENT, YCENT;
  float activ;
  int ID;
  float divide;
  int num;
  float tx = 0, ty = 0, tz = 0;
  float xt = 0, yt = 0, zt = 0;
  float xr = 0, yr = 0, zr = 0;
  float prev_motion[6];
  float motion_phase_index = 0.0, motion_incr_ph_index = 0.0;
  float time_per_frame = 0;

  POINT C;

  for(i = 0; i < 6; i++)
    prev_motion[i] = 0.0;
/*
** -------------------
** READ COMMAND LINE :
** -------------------
*/
  while (--argc > 0 && (*++argv)[0] == '-') {
    char    *s;
    for (s = argv[0] + 1; *s; s++)
      switch (*s) {
        case 'x':
          xt = atof(*++argv);
          motion_flag = 1;
          argc--;
          break;
        case 'y':
          yt = atof(*++argv);
          motion_flag = 1;
          argc--;
          break;
        case 'z':
          zt = atof(*++argv);
          motion_flag = 1;
          argc--;
          break;
        case 'X':
          xr = atof(*++argv);
          motion_flag = 1;
          argc--;
          break;
        case 'Y':
          yr = atof(*++argv);
          motion_flag = 1;
          argc--;
          break;
        case 'Z':
          zr = atof(*++argv);
          motion_flag = 1;
          argc--;
          break;
        default:
          fprintf (stderr, "Unknown option: -%c\n", *s);
          Abort ("");
          break;
     }
   }

  if (argc < 2) {
    fprintf(stderr, DOCUMENTATION, argv[0]);
    exit(1);
  }
  else {
    if (sscanf(*argv, "%s",gen_parfile)==0)
      Abort("Can not get gen_parfile filename");
    if(sscanf(*++argv, "%s", outputbase)==0)          
      Abort("Can not get output base name");
  }

/*
** -----------------
** OPEN *_LOG FILE :
** -----------------
*/
  strcpy(logfile, outputbase);
  strcat(logfile, "_log");
  if((fp = fopen(logfile, "w")) == NULL)
    Abort("Can not open *_log output file");
  fprintf(fp,"\n*** DYNAMIC NCAT BRAIN PHANTOM, version 1.0 in C***");
  fprintf(fp,"\n--------------------------------------------");

/*
** -------------------
** READ PARAMETER FILE
** -------------------
*/
  GET_DYN_PARAMS(gen_parfile);
  if(out_frames > 1)
    time_per_frame = out_period / (out_frames-1);
  else
    time_per_frame = 0;

  if(motion_flag)
    {
    x_trans = xt;
    y_trans = yt;
    z_trans = zt;
    x_rot = xr;
    y_rot = yr;
    z_rot = zr;
    }


/* set phantom dimensions: */
  XDIM=array_size;
  YDIM=array_size;
  ZDIM_OUTPUT=1+(endslice-startslice);

/*
** ---------------------------
** DYNAMICALLY ALLOCATE MEMORY
** ---------------------------
**
** arrays holding output phantoms:
*/
  txdim = subvoxel_index * XDIM + 10 * subvoxel_index;
  tydim = subvoxel_index * YDIM + 10 * subvoxel_index;
  tzdim = ZDIM_OUTPUT * subvoxel_index;
  
  xoffset = 10 * subvoxel_index;
  yoffset = 10 * subvoxel_index;
  zoffset = 0;

  fphan =    vector(0, XDIM*YDIM*ZDIM_OUTPUT-1);
  fphan_id =    vector(0, XDIM*YDIM*ZDIM_OUTPUT-1);
  fphan_avg =    vector(0, XDIM*YDIM*ZDIM_OUTPUT-1);

  slice = vector(0, txdim*tydim);

/*
** OFFSETS: Used to center torso in the final 3D image produced
*/
  XCENT=(XDIM/2.0);
  YCENT=(YDIM/2.0);

/*
** VOXEL VOLUME in mL (cm^3)
*/
  voxel_volume=pixel_width * pixel_width * slice_width;
  subvxl_vol=(pixel_width * pixel_width * slice_width) / (subvoxel_index*subvoxel_index*subvoxel_index);

/*
** ----------------------
** READ TORSO NURB FILE
** ----------------------
*/
  xoff = 0.0; yoff = 0.0; zoff = 0.0;
  centerx = 93.0; centery = 106.0;  //as determined from skull in Rhino3d 
  centerx /= (pixel_width*10.0);
  centery /= (pixel_width*10.0);

  xoff = (XCENT - centerx) * pixel_width * 10.0;
  yoff = (YCENT - centery) * pixel_width * 10.0;
  zoff= -( (startslice-1)*slice_width*10.0) + 10;

 
  for(i = 0; i <= END_MODELS; i++)
  {
    tmodel[i].num_tris = 0;
    sprintf(filename, "vtk/t%i.vtk", i);
    Read_Triangle_ModelVTK(filename, &tmodel[i]);
  }

  slice_tmodel.tris = tri_vector(0, 20000);

  /* Read motion curve */
  if(user_time_flag == 0)
    {
    printf("\nUsing default motion curve");
    for(i = 0; i < 6; i++)
      Read_motion_file("brain_curve.nrb", &motion_curve[i]);
    }

  motion_incr_ph_index = time_per_frame/motion_period;


/*
for(i = 0; i < 100; i++)
  {
  printf("\n");
  for(j = 0; j < 6; j++)
    {
    C = CurvePoint(motion_curve[j].pol.n-1, 3, motion_curve[j].knt, motion_curve[j].pol.Pw, (float)i/100.0);
    printf(" %f", C.y);
    }
  }
exit(1);
*/


/*
** Initialize the arrays containing the time average phantoms:
*/
   for(i = 0; i <= XDIM*YDIM*ZDIM_OUTPUT-1; i++)
     fphan_avg[i] =0.;


  tx = (tmodel[21].max_x + tmodel[21].min_x)/2.0; //approximate neck point
  ty = (tmodel[21].max_y + tmodel[21].min_y)/2.0;
  tz = tmodel[21].min_z * (slice_width/subvoxel_index*10.0);

  if(user_time_flag == 1)
    {
    x_trans = 1;
    y_trans = 1;
    z_trans = 1;
    x_rot = 1;
    y_rot = 1;
    z_rot = 1;
    }

  if(motion_option == 0)
    {
    C = CurvePoint(motion_curve[0].pol.n-1, 3, motion_curve[0].knt, motion_curve[0].pol.Pw, motion_start_ph_index);
    xr = C.y * x_rot;

    C = CurvePoint(motion_curve[1].pol.n-1, 3, motion_curve[1].knt, motion_curve[1].pol.Pw, motion_start_ph_index);
    yr = C.y * y_rot;

    C = CurvePoint(motion_curve[2].pol.n-1, 3, motion_curve[2].knt, motion_curve[2].pol.Pw, motion_start_ph_index);
    zr = C.y * z_rot;

    C = CurvePoint(motion_curve[3].pol.n-1, 3, motion_curve[3].knt, motion_curve[3].pol.Pw, motion_start_ph_index);
    xt = C.y * x_trans;

    C = CurvePoint(motion_curve[4].pol.n-1, 3, motion_curve[4].knt, motion_curve[4].pol.Pw, motion_start_ph_index);
    yt = C.y * y_trans;

    C = CurvePoint(motion_curve[5].pol.n-1, 3, motion_curve[5].knt, motion_curve[5].pol.Pw, motion_start_ph_index);
    zt = C.y * z_trans;

    for(i = 0; i <= END_MODELS; i++)
      MoveTriangleModel(&tmodel[i], xt, yt, zt, xr, yr, zr, &tx, &ty, &tz);
    }

  for(kf=1; kf <= out_frames; kf++) {
    fprintf(fp,"\5--------------------------------------------------");
    fprintf(fp,"\n--------------------------------------------------");
    fprintf(fp,"\nCREATING FRAME #%i ...\n",kf);
    printf("\nCREATING FRAME #%i ...\n",kf);
    fprintf(fp,"\nTime = %7.3f seconds",time_per_frame*(kf-1));
 
    if(motion_option == 1)
      {
      motion_phase_index = motion_start_ph_index + motion_incr_ph_index*(kf-1);
      while (motion_phase_index>=1.000)
        motion_phase_index = motion_phase_index - 1;

      C = CurvePoint(motion_curve[0].pol.n-1, 3, motion_curve[0].knt, motion_curve[0].pol.Pw, motion_phase_index);
      xr = (C.y-prev_motion[0]) * x_rot;
      prev_motion[0] = C.y;
      printf("\nMOTION for frame %i, time = %f: angx = %f ", kf, C.x, C.y*x_rot);

      C = CurvePoint(motion_curve[1].pol.n-1, 3, motion_curve[1].knt, motion_curve[1].pol.Pw, motion_phase_index);
      yr = (C.y-prev_motion[1]) * y_rot;
      prev_motion[1] = C.y;
      printf("angy = %f ", C.y*y_rot);

      C = CurvePoint(motion_curve[2].pol.n-1, 3, motion_curve[2].knt, motion_curve[2].pol.Pw, motion_phase_index);
      zr = (C.y-prev_motion[2]) * z_rot;
      prev_motion[2] = C.y;
      printf("angz = %f ", C.y*z_rot);

      C = CurvePoint(motion_curve[3].pol.n-1, 3, motion_curve[3].knt, motion_curve[3].pol.Pw, motion_phase_index);
      xt = (C.y-prev_motion[3]) * x_trans;
      prev_motion[3] = C.y;
      printf("tx = %f ", C.y*x_trans);

      C = CurvePoint(motion_curve[4].pol.n-1, 3, motion_curve[4].knt, motion_curve[4].pol.Pw, motion_phase_index);
      yt = (C.y-prev_motion[4]) * y_trans;
      prev_motion[4] = C.y;
      printf("ty = %f ", C.y*y_trans);

      C = CurvePoint(motion_curve[5].pol.n-1, 3, motion_curve[5].knt, motion_curve[5].pol.Pw, motion_phase_index);
      zt = (C.y-prev_motion[5]) * z_trans;
      prev_motion[5] = C.y;
      printf("tz = %f", C.y*z_trans);

      for(i = 0; i <= END_MODELS; i++)
        MoveTriangleModel(&tmodel[i], xt, yt, zt, xr, yr, zr, &tx, &ty, &tz);
      }
    else
      motion_phase_index = motion_start_ph_index;

    fprintf(fp,"\nCurrent motion phase index = %7.3f ",motion_phase_index);
    fprintf(fp,"\nTranslating by %f %f %f", xt, yt, zt);
    fprintf(fp,"\nRotating by %f %f %f", xr, yr, zr);

/*
**  Initialize time frame arrays containing the phantoms:
*/
   for(i = 0; i <= XDIM*YDIM*ZDIM_OUTPUT-1; i++)
     {
     fphan[i] =0.;
     fphan_id[i] = 0.;
     }

//sprintf(filename, "body%i.raw", kf);
//Output_Triangle_Model(filename, &tmodel[21]);

   printf("\nRendering....\n");
/*
**  -------------------
**  | RENDER LOOP     |
**  -------------------
*/
   for(j = 0; j < tzdim; j++) {
     printf("\rSlice = %d/%d (%.1f%%)", j, tzdim, j * 100.0 / tzdim);
     /*Initialize slice image */
     for(i = 0; i < txdim*tydim; i++)
       slice[i] = -1.0;

     //Render Scalp first
     k = 21;
     intensity = k;
     num = 0;
     slice_tmodel.num_tris = num;
     for(i = 0; i < tmodel[k].num_tris; i++)
	 {
       if(Test_extents_tri_z(tmodel[k].tris[i], j))
       {
         slice_tmodel.tris[num].vertex[0].x = tmodel[k].tris[i].vertex[0].x;
         slice_tmodel.tris[num].vertex[0].y = tmodel[k].tris[i].vertex[0].y;
         slice_tmodel.tris[num].vertex[0].z = tmodel[k].tris[i].vertex[0].z;

         slice_tmodel.tris[num].vertex[1].x = tmodel[k].tris[i].vertex[1].x;
         slice_tmodel.tris[num].vertex[1].y = tmodel[k].tris[i].vertex[1].y;
         slice_tmodel.tris[num].vertex[1].z = tmodel[k].tris[i].vertex[1].z;

         slice_tmodel.tris[num].vertex[2].x = tmodel[k].tris[i].vertex[2].x;
         slice_tmodel.tris[num].vertex[2].y = tmodel[k].tris[i].vertex[2].y;
         slice_tmodel.tris[num].vertex[2].z = tmodel[k].tris[i].vertex[2].z;

         num++;
         slice_tmodel.num_tris = num;
       }
	 }
     if(slice_tmodel.num_tris > 0)
     {
       for(i = 0; i < slice_tmodel.num_tris; i++)
         Render_triangle(slice_tmodel.tris[i], intensity, slice, j);
       Fill_tri(slice_tmodel, intensity, slice, j);
     }

     //Render Skull second
     k = 28;
     intensity = k;
     num = 0;
     slice_tmodel.num_tris = num;
     for(i = 0; i < tmodel[k].num_tris; i++)
	 {
       if(Test_extents_tri_z(tmodel[k].tris[i], j))
       {
         slice_tmodel.tris[num].vertex[0].x = tmodel[k].tris[i].vertex[0].x;
         slice_tmodel.tris[num].vertex[0].y = tmodel[k].tris[i].vertex[0].y;
         slice_tmodel.tris[num].vertex[0].z = tmodel[k].tris[i].vertex[0].z;

         slice_tmodel.tris[num].vertex[1].x = tmodel[k].tris[i].vertex[1].x;
         slice_tmodel.tris[num].vertex[1].y = tmodel[k].tris[i].vertex[1].y;
         slice_tmodel.tris[num].vertex[1].z = tmodel[k].tris[i].vertex[1].z;

         slice_tmodel.tris[num].vertex[2].x = tmodel[k].tris[i].vertex[2].x;
         slice_tmodel.tris[num].vertex[2].y = tmodel[k].tris[i].vertex[2].y;
         slice_tmodel.tris[num].vertex[2].z = tmodel[k].tris[i].vertex[2].z;

         num++;
         slice_tmodel.num_tris = num;
       }
	 }
     if(slice_tmodel.num_tris > 0)
     {
       for(i = 0; i < slice_tmodel.num_tris; i++)
         Render_triangle(slice_tmodel.tris[i], intensity, slice, j);
       Fill_tri(slice_tmodel, intensity, slice, j);
     }

     //Render Spinal fluid third
     k = 255;
     intensity = k;
     num = 0;
     slice_tmodel.num_tris = num;
     for(i = 0; i < tmodel[k].num_tris; i++)
	 {
       if(Test_extents_tri_z(tmodel[k].tris[i], j))
       {
         slice_tmodel.tris[num].vertex[0].x = tmodel[k].tris[i].vertex[0].x;
         slice_tmodel.tris[num].vertex[0].y = tmodel[k].tris[i].vertex[0].y;
         slice_tmodel.tris[num].vertex[0].z = tmodel[k].tris[i].vertex[0].z;

         slice_tmodel.tris[num].vertex[1].x = tmodel[k].tris[i].vertex[1].x;
         slice_tmodel.tris[num].vertex[1].y = tmodel[k].tris[i].vertex[1].y;
         slice_tmodel.tris[num].vertex[1].z = tmodel[k].tris[i].vertex[1].z;

         slice_tmodel.tris[num].vertex[2].x = tmodel[k].tris[i].vertex[2].x;
         slice_tmodel.tris[num].vertex[2].y = tmodel[k].tris[i].vertex[2].y;
         slice_tmodel.tris[num].vertex[2].z = tmodel[k].tris[i].vertex[2].z;

         num++;
         slice_tmodel.num_tris = num;
       }
	 }
	 if(slice_tmodel.num_tris > 0)
     {
       for(i = 0; i < slice_tmodel.num_tris; i++)
         Render_triangle(slice_tmodel.tris[i], intensity, slice, j);
       Fill_tri(slice_tmodel, intensity, slice, j);
     }

     for(k = 0; k <= END_MODELS; k++)
     {
       if(k != 21 && k != 28 && k != 255)
       {
         intensity = k;
         num = 0;
         slice_tmodel.num_tris = num;       
         for(i = 0; i < tmodel[k].num_tris; i++)
         {
           if(Test_extents_tri_z(tmodel[k].tris[i], j))
           {
             slice_tmodel.tris[num].vertex[0].x = tmodel[k].tris[i].vertex[0].x;
             slice_tmodel.tris[num].vertex[0].y = tmodel[k].tris[i].vertex[0].y;
             slice_tmodel.tris[num].vertex[0].z = tmodel[k].tris[i].vertex[0].z;

             slice_tmodel.tris[num].vertex[1].x = tmodel[k].tris[i].vertex[1].x;
             slice_tmodel.tris[num].vertex[1].y = tmodel[k].tris[i].vertex[1].y;
             slice_tmodel.tris[num].vertex[1].z = tmodel[k].tris[i].vertex[1].z;

             slice_tmodel.tris[num].vertex[2].x = tmodel[k].tris[i].vertex[2].x;
             slice_tmodel.tris[num].vertex[2].y = tmodel[k].tris[i].vertex[2].y;
             slice_tmodel.tris[num].vertex[2].z = tmodel[k].tris[i].vertex[2].z;

             num++;
             slice_tmodel.num_tris = num;
           }
         }
         if(slice_tmodel.num_tris > 0)
         {
           for(i = 0; i < slice_tmodel.num_tris; i++)
             Render_triangle(slice_tmodel.tris[i], intensity, slice, j);
           Fill_tri(slice_tmodel, intensity, slice, j);
         }
       }
     }

     for(k = xoffset; k < txdim; k++)
       for(l = yoffset; l < tydim; l++)
       {
         index1 = (int)( (k-xoffset) / subvoxel_index) + 
                  (int)( (l-yoffset) / subvoxel_index) * XDIM + (int)(j/subvoxel_index) * XDIM * YDIM;
         index2 = k + l * txdim;
        
         if(slice[index2] != -1)
         {
           ID = slice[index2];
           activ = activity[ID];
           fphan_id[index1] += slice[index2] / (subvoxel_index*subvoxel_index*subvoxel_index);
           fphan[index1] += activ / (subvoxel_index*subvoxel_index*subvoxel_index);
           fphan_avg[index1] += activ / (subvoxel_index*subvoxel_index*subvoxel_index*out_frames);
         }
       }
    } /*End slice loop*/
    printf("\rSlice = %d/%d (100%%) \n", tzdim, tzdim);

//  GetVolumes(fphan_id); 
    if (act_phan_each == 1) 
    { /* SAVE ACTIVITY PHANTOM for frame kf to file */
      /* write activity phantom for frame kf */
      sprintf(outfile,"%s_%d",outputbase,kf);
      fprintf(fp, "\nSaving activity phantom in file : %s", outfile);
      SAVE_TO_FILE(fphan,XDIM,YDIM,ZDIM_OUTPUT,outfile);
    }
  } /* End kf loop*/

/*
    --------------------------
    | END of MAIN VOXEL LOOP |
    --------------------------
*/
  if (act_phan_ave == 1) 
  {/* SAVE AVERAGE ACTIVITY PHANTOM to file */
    sprintf(outfile,"%s_%s",outputbase,"av");
    fprintf(fp, "\nSaving time-avg activity phantom in file : %s", outfile);
    fprintf(fp, "\n-------------------------------------------------");
    SAVE_TO_FILE(fphan_avg,XDIM,YDIM,ZDIM_OUTPUT,outfile);    
  }
  fclose(fp);
} /*	End Main	*/
