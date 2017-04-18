#include "global_includes.h"
#include "dncatsubs.h"
#include "nurbs.h"
#include "global_vars.h"
#include "constants.h"

#include <algorithm>

/*--------------------------------------------------------------*/   
/*--------------------------------------------------------------*/   
void NEXTLINE(FILE *fp)
/*----------------------------------------
** advance file pointer (fp) to next line.
**----------------------------------------
*/
{    
 	int c = 0;
    	while(c != '\n' && c != EOF)
          c = fgetc(fp);
}   /* end NEXTLINE subroutine*/ 



float Calc_Angle(float dx, float dy)
{
  float angle = 0.0;

  if(dx > 0 && dy > 0)
    angle = 180.0/PI*atan(dy/dx);
  else if(dx < 0 && dy > 0)
    angle = 180.0 + 180.0/PI*atan(dy/dx);
  else if(dx < 0 && dy < 0)
    angle = 180.0 + 180.0/PI*atan(dy/dx);
  else if(dx > 0 && dy < 0)
    angle = 360.0 + 180.0/PI*atan(dy/dx);                                      
  else if(dx == 0 && dy == 0)
    angle = 0.0;
  else if(dx == 0 && dy > 0)
    angle = 90.0;
  else if(dx == 0 && dy < 0)
    angle = 270.0;
  else if(dx > 0 && dy == 0)
    angle = 0.0;
  else if(dx < 0 && dy == 0)
    angle = 180.0;
  
  return angle;
}

/*--------------------------------------------------------------*/   
/*--------------------------------------------------------------*/   
float sqr(float x)
/*---------------------------------------------------------------
**  This subroutine performs the square function
**----------------------------------------------------------------
*/
{
  return(x*x);
}
/*--------------------------------------------------------------*/   


void READ_BRAIN_CURVES(char *filename, CURVE motion_curve[6])
{
  FILE  *fp;
  int num, flag, count;
  // char comma;
  float time, ang_x, ang_y, ang_z, t_x, t_y, t_z, period;
  int i;

/*
** READ PARAMETER FILE
*/
  if((fp = fopen(filename, "r")) == NULL )
    Abort("Cannot open brain curves file");

  /* Read file to determine how many points there are in total */
  NEXTLINE(fp);
  count = 0;
 flag = fscanf(fp, "%f %f %f %f %f %f %f", &time, &ang_x, &ang_y, &ang_z, &t_x, &t_y, &t_z);

  while(flag)
    {
    count++;
    period = time;
 flag = fscanf(fp, "%f %f %f %f %f %f %f", &time, &ang_x, &ang_y, &ang_z, &t_x, &t_y, &t_z);
    }
  num = count;

  printf("\n  Total of %i points determine the brain motion curves", num);
  printf("\n  Duration of curve = %f seconds", period);
  motion_period = period;

  for(i = 0; i < 6; i++)
    {
    MakeCurve(&motion_curve[i], num, num+p_degree, p_degree);
    Calc_UniformKnotVector(motion_curve[i].pol.n-1, p_degree, motion_curve[i].knt, motion_curve[i].uk);
    }
  fclose(fp);

  if((fp = fopen(filename, "r")) == NULL )
    Abort("Cannot open brain curves file");

  NEXTLINE(fp);
  count = 0; 
  flag = fscanf(fp, "%f %f %f %f %f %f %f", &time, &ang_x, &ang_y, &ang_z, &t_x, &t_y, &t_z);
  while(flag)
    {
    motion_curve[0].pol.Pw[count].x = time;
    motion_curve[0].pol.Pw[count].y = ang_x;
    motion_curve[0].pol.Pw[count].z = period;
    
    motion_curve[1].pol.Pw[count].x = time;
    motion_curve[1].pol.Pw[count].y = ang_y;
    motion_curve[1].pol.Pw[count].z = period;

    motion_curve[2].pol.Pw[count].x = time;
    motion_curve[2].pol.Pw[count].y = ang_z;
    motion_curve[2].pol.Pw[count].z = period;

    motion_curve[3].pol.Pw[count].x = time;
    motion_curve[3].pol.Pw[count].y = t_x;
    motion_curve[3].pol.Pw[count].z = period;

    motion_curve[4].pol.Pw[count].x = time;
    motion_curve[4].pol.Pw[count].y = t_y;
    motion_curve[4].pol.Pw[count].z = period;

    motion_curve[5].pol.Pw[count].x = time;
    motion_curve[5].pol.Pw[count].y = t_z;
    motion_curve[5].pol.Pw[count].z = period;

    count++;
    flag = fscanf(fp, "%f %f %f %f %f %f %f", &time, &ang_x, &ang_y, &ang_z, &t_x, &t_y, &t_z);
    }

  if (num < 3) {
    printf("\n!!! Number of time points for curves must be greater than 3 for cubic interpolation");
    exit(1);
  }
  
  fclose(fp);
}

/*--------------------------------------------------------------*/   
/*--------------------------------------------------------------*/   
void GET_DYN_PARAMS(char *parfile)
/*---------------------------------------------------------------
**  This subroutine reads the parameter file and checks the
**  parameter values for validity. ALL params are global variables.
**  If an invalid parameter value is discovered, the bad_param flag
**  is set to 1 and program is terminated upon return to main part.
**----------------------------------------------------------------
*/
{
  FILE  *fp;

/* 
** READ PARAMETER FILE
*/
  if((fp = fopen(parfile, "r")) == NULL )
    Abort("Cannot open parfile");

  fscanf(fp, "%i", &act_phan_each);             NEXTLINE(fp);
  fscanf(fp, "%i", &act_phan_ave);              NEXTLINE(fp);

  fscanf(fp, "%i", &motion_option);             NEXTLINE(fp);

  fscanf(fp, "%f", &out_period);                NEXTLINE(fp);
  fscanf(fp, "%i", &out_frames);                NEXTLINE(fp);

  fscanf(fp, "%f", &motion_period);             NEXTLINE(fp);
  fscanf(fp, "%f", &motion_start_ph_index);     NEXTLINE(fp);
  fscanf(fp, "%i", &user_time_flag);		NEXTLINE(fp);
  fscanf(fp, "%f", &x_trans);			NEXTLINE(fp);
  fscanf(fp, "%f", &y_trans);                   NEXTLINE(fp);
  fscanf(fp, "%f", &z_trans);                   NEXTLINE(fp);
  fscanf(fp, "%f", &x_rot);                     NEXTLINE(fp);
  fscanf(fp, "%f", &y_rot);                     NEXTLINE(fp);
  fscanf(fp, "%f", &z_rot);                     NEXTLINE(fp);
  fscanf(fp, "%s", user_motion_file);		NEXTLINE(fp);
  fscanf(fp, "%f", &pixel_width);  		NEXTLINE(fp);
  fscanf(fp, "%f", &slice_width);  		NEXTLINE(fp);
  fscanf(fp, "%i", &array_size);  		NEXTLINE(fp);
  fscanf(fp, "%i", &startslice); 		NEXTLINE(fp);
  fscanf(fp, "%i", &endslice); 			NEXTLINE(fp);
  fscanf(fp, "%i", &subvoxel_index);            NEXTLINE(fp);

fscanf(fp, "%f", &activity[10]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[70]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[2]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[50]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[75]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[15]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[5]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[80]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[6]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[90]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[1]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[85]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[7]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[27]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[114]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[9]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[88]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[52]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[60]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[41]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[19]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[159]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[32]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[56]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[110]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[74]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[145]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[61]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[130]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[64]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[140]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[164]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[26]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[62]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[165]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[119]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[99]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[196]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[139]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[118]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[125]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[18]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[132]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[251]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[38]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[98]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[63]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[154]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[97]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[37]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[175]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[54]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[112]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[69]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[4]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[108]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[53]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[39]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[16]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[14]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[11]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[12]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[203]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[102]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[133]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[255]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[232]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[233]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[8]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[3]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[20]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[13]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[17]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[30]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[45]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[73]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[105]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[57]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[59]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[83]); NEXTLINE(fp);   
fscanf(fp, "%f", &activity[243]); NEXTLINE(fp);
fscanf(fp, "%f", &activity[21]); NEXTLINE(fp);
fscanf(fp, "%f", &activity[28]); NEXTLINE(fp);
fscanf(fp, "%f", &activity[24]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[36]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[101]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[25]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[72]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[254]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[29]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[35]); NEXTLINE(fp);     
fscanf(fp, "%f", &activity[34]); NEXTLINE(fp);     
fscanf(fp, "%f", &activity[23]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[33]); NEXTLINE(fp);  
fscanf(fp, "%f", &activity[128]); NEXTLINE(fp);     
fscanf(fp, "%f", &activity[43]); NEXTLINE(fp);     
fscanf(fp, "%f", &activity[76]); NEXTLINE(fp); 
fscanf(fp, "%f", &activity[67]); NEXTLINE(fp); 


  if(user_time_flag == 1)
    {
    printf("\nBrain motion determined by a user defined time curve");
    READ_BRAIN_CURVES(user_motion_file, motion_curve);
    }

  fclose(fp);
} /* end get_params subroutine */
/*--------------------------------------------------------------*/   

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
void Read_motion_file(char *file, CURVE *C)
/*------------------------------------------------------------------
**  This subroutine is used to read in motion curves for the
**  respiratory movements
**------------------------------------------------------------------
*/
{
  FILE *m_fp;
  int i, n;
  int p = 3;
  char line[30];
  float temp, tmp_x, tmp_y, tmp_z;

  if((m_fp = fopen(file, "r")) == NULL)
    Abort("Can not open motion file");

  /*Read in number of curve points*/
  fscanf(m_fp, "%i", &n);   fscanf(m_fp, "%s", line);

  /*Create memory for curve*/
  MakeCurve(C, n, n+p, p);

  /*Read in Knot Vector*/
  fscanf(m_fp, "%s", line); fscanf(m_fp, "%s", line);
  for(i = 0; i <= n+p; i++)
    {
    fscanf(m_fp, "%f", &temp);
    C->knt.U[i] = temp;
    }

  /*Read in Control Points*/
  fscanf(m_fp, "%s", line); fscanf(m_fp, "%s", line);
  for(i = 0; i < n; i++)
    {
    fscanf(m_fp, "%f %f %f", &tmp_x, &tmp_y, &tmp_z);
    C->pol.Pw[i].x = tmp_x;
    C->pol.Pw[i].y = tmp_y;
    C->pol.Pw[i].z = tmp_z;
    C->pol.Pw[i].w = 0.0;
    }

  fclose(m_fp);
}

/*------------------------------------------------------------------*/   
/*------------------------------------------------------------------*/   
void vm_mult(float m[4][4], float p[4], float result[4])
/*-------------------------------------------------------------------
**  This subroutine performs the multiplication of a vector (p) with
**  a matrix (m) and returns the result (result).
**-------------------------------------------------------------------
*/
{
  result[1] = m[1][1] * p[1] + m[2][1] * p[2] + m[3][1] * p[3];
  result[2] = m[1][2] * p[1] + m[2][2] * p[2] + m[3][2] * p[3];
  result[3] = m[1][3] * p[1] + m[2][3] * p[2] + m[3][3] * p[3];
}
/*-------------------------------------------------------------------*/   


/*---------------------------------------------------------------------------------------*/   
/*---------------------------------------------------------------------------------------*/   
void Rotate(float x_rot, float y_rot, float z_rot, float *p, float tx, float ty, float tz)
/*---------------------------------------------------------------------------------------
**  This subroutine is used to rotate a point (p) about the x, y 
**  or z axes.
**  tx, ty, tz = pivot point  
**  x_rot = angle to rotate about the x-axis
**  y_rot = angle to rotate about the y-axis
**  z_rot = angle to rotate about the z-axis
**---------------------------------------------------------------------------------------
*/
{
  float result[4];
  float Rx[4][4], Rz[4][4], Ry[4][4];
  int i, j;
  float cosine, sine;

  for(i = 1; i<=3; i++)
    for(j = 1; j<=3; j++)
      {
      Rx[i][j] = 0.0;
      Rz[i][j] = 0.0;
      Ry[i][j] = 0.0;
      }
    
  cosine = cos(x_rot * PI / 180);
  sine = sin(x_rot * PI / 180);
    
  Rx[1][1] = 1.0;
  Rx[2][2] = cosine;
  Rx[2][3] = -sine;
  Rx[3][2] = sine;
  Rx[3][3] = cosine;
  
  cosine = cos(z_rot * PI / 180);
  sine = sin(z_rot * PI / 180);

  Rz[1][1] = cosine; 
  Rz[1][2] = -sine;
  Rz[2][1] = sine;
  Rz[2][2] = cosine;
  Rz[3][3] = 1.0;

  cosine = cos(y_rot * PI / 180);
  sine = sin(y_rot * PI / 180);
  
  Ry[1][1] = cosine;
  Ry[1][3] = sine;   
  Ry[2][2] = 1.0;  
  Ry[3][1] = -sine;
  Ry[3][3] = cosine;   

/* Translate pivot to origin */
  p[1] = p[1] - tx;      
  p[2] = p[2] - ty;      
  p[3] = p[3] - tz;
      
/* Rotate about x */
  vm_mult(Rx, p, result);
  p[1] = result[1]; p[2] = result[2]; p[3] = result[3];

/* Rotate about y */
  vm_mult(Ry, p, result);
  p[1] = result[1]; p[2] = result[2]; p[3] = result[3];

/* Rotate about z */
  vm_mult(Rz, p, result);
  p[1] = result[1]; p[2] = result[2]; p[3] = result[3];


/* Translate back */
  p[1] = p[1] + tx;      
  p[2] = p[2] + ty;      
  p[3] = p[3] + tz;
}
/*---------------------------------------------------------------------------------------*/   

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
void Output2(SURFACE *nrb_model)
/*------------------------------------------------------------------
**  This subroutine is used to output a NURBS surface to a file.
**  This was used to debug the program.
**------------------------------------------------------------------
*/
{
  int i, j;
  float pt[4];

  printf("\nSrfControlPtGrid");
  printf("\n3\n3");
  printf("\n%i",nrb_model->net.m);
  printf("\n%i",nrb_model->net.n);

  for(j = 0; j < nrb_model->net.n; j++)
    {
    for(i = 0; i < nrb_model->net.m; i++)
      {
      pt[1] = nrb_model->net.Pw[j][i].x;
      pt[2] = nrb_model->net.Pw[j][i].y;
      pt[3] = nrb_model->net.Pw[j][i].z;

    printf("\n%f,%f,%f",pt[1],pt[2],pt[3]);
    }
  }
}
/*----------------------------------------------------------------*/

/*----------------------------------------------------------------*/
/*----------------------------------------------------------------*/
float distance(float x1, float y1, float z1, float x2, float y2, float z2)
/*----------------------------------------------------------------------
**  This subroutine is used to compute the distance between two points.
**----------------------------------------------------------------------
*/
{
float result;

x1 = x1 - x2;
y1 = y1 - y2;
z1 = z1 - z2;

result = sqrt(x1*x1 + y1*y1 + z1*z1);
return result;
}


/*-----------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------*/
void scale(float x_center, float y_center, float z_center, float xbegin, float ybegin,
float zbegin, float *xend, float *yend, float *zend, float factor)
/*------------------------------------------------------------------------------------
**  This subroutine is used to scale a point (xbegin, ybegin, zbegin) by the amount
**  factor.  The origin for the scaling operation is the point 
**  (x_center, y_center, z_center). The scaled point is returned as (xend, yend, zend).
**------------------------------------------------------------------------------------
*/
{
  *xend = xbegin;
  *yend = ybegin;
  *zend = zbegin;

  *xend -= x_center;
  *yend -= y_center;
  *zend -= z_center;

  *xend *= factor;
  *yend *= factor;
  *zend *= factor;

  *xend += x_center;
  *yend += y_center;
  *zend += z_center;
}

/*--------------------------------------------------------------------------------------------*/   
float Set_scale(float xp, float yp, float zp, float *xend, float *yend, float *zend, float xc,
float yc, float zc, float dist, float x1, float x2, float xacc)
/*--------------------------------------------------------------------------------------------*/   
/*--------------------------------------------------------------------------------------------
**  This subroutine sets the distance between two points by scaling
**  in 2D
**--------------------------------------------------------------------------------------------
*/
{
  int j;
  float dx, f, fmid, xmid, rtb;
  float tmp_dist;

  scale(xc, yc, zc, xp, yp, zp, xend, yend, zend, x1);
  tmp_dist= distance(xp, yp, zp, *xend, *yend, *zend);

  f = dist - tmp_dist;
  if(fabs(f) < xacc)
    return x1;

  scale(xc, yc, zc, xp, yp, zp, xend, yend, zend, x2);
  tmp_dist= distance(xp, yp, zp, *xend, *yend, *zend);

  fmid = dist - tmp_dist;
  if(fabs(fmid) < xacc)
    return x2;

  if(f*fmid >= 0.0)
    {
    return 0;
    }

  rtb = f < 0.0 ? (dx= x2 - x1, x1) : (dx= x1-x2, x2);

  while(fabs(fmid) > xacc && fmid != 0.0)
    {
    dx *= 0.5;
    xmid = rtb + dx;

    scale(xc, yc, zc, xp, yp, zp, xend, yend, zend, xmid);
    tmp_dist= distance(xp, yp, zp, *xend, *yend, *zend);

    fmid = dist - tmp_dist;
    
    if(fmid <= 0.0)
      rtb = xmid;  
    }
  return rtb; 
}
   

/*----------------------------------------------------------------*/   
void Calc_extents(SURFACE *nrb_model)
/*----------------------------------------------------------------*/   
/*------------------------------------------------------------------
**  This subroutine is used to set the maximum and minimum dimensions 
**  of a given NURBS model.  These dimensions are used to pre-render
**  the NURBS surface into its own image before combining the image
**  into the final output image.
**------------------------------------------------------------------
*/
{
  int i, j;
  float tz;
  float max_z, min_z;

  /*Initialize the maximum and minimum to the first control point*/
  max_z = nrb_model->net.Pw[0][0].z;
  min_z = nrb_model->net.Pw[0][0].z;

  /*Cycle through the model's control points*/
  for(i = 0; i < nrb_model->net.m; i++)  {
    for(j = 0; j< nrb_model->net.n; j++)  {
      tz = nrb_model->net.Pw[j][i].z;

      /*Find the maximum*/
      if(tz > max_z)
        max_z = tz;

      /*Find the minimum*/
      if(tz < min_z)
        min_z = tz;
    }
  }
  nrb_model->max_z = floor( (max_z) / (slice_width/subvoxel_index * 10.0) + 0.5);
  nrb_model->min_z = floor( (min_z) / (slice_width/subvoxel_index * 10.0) + 0.5);
}
/*----------------------------------------------------------------*/   

/*--------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------*/
void initialize_structure(SURFACE *nrb_model, int n, int m)
/*------------------------------------------------------------------------------
**  This subroutine initializes the specified NURBS surface for the heart
**    id:       Heart surface to be initialized
**    m,n:      number of control points in the v and u directions of the surface
**-------------------------------------------------------------------------------
*/
{
  int i, j;

  /*Setup NURB surfaces of the heart*/
  nrb_model->net.Pw = cp_matrix(0, n-1, 0 , m-1);
  nrb_model->net.n = n;
  nrb_model->net.m = m;
  nrb_model->knu.U = vector(0, n + p_degree);
  nrb_model->knu.m = n + p_degree;
  nrb_model->knv.U = vector(0, m + q_degree);
  nrb_model->knv.m = m + q_degree;

  nrb_model->flag = 0;

  for(i = 0; i < m; i++)
    for(j = 0; j < n; j++)
      {
      nrb_model->net.Pw[j][i].x = 0.0;
      nrb_model->net.Pw[j][i].y = 0.0;
      nrb_model->net.Pw[j][i].z = 0.0;
      nrb_model->net.Pw[j][i].w = 0.0;
      }
}
/*-------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------*/   
/*-------------------------------------------------------------------------------*/   
void initialize_structure2(SURFACE *nrb_model, int n, int m)
/*------------------------------------------------------------------------------
**  This subroutine initializes the specified NURBS surface for the heart
**    id:	Heart surface to be initialized
**    m,n:	number of control points in the v and u directions of the surface 
**-------------------------------------------------------------------------------
*/
{
  int i, j;

  /*Setup NURB surfaces of the heart*/
  nrb_model->net.Pw = cp_matrix(0, n-1, 0 , m-1);
  nrb_model->net.n = n;
  nrb_model->net.m = m;
  nrb_model->knu.U = vector(0, n + p_degree);
  nrb_model->knu.m = n + p_degree;
  nrb_model->knv.U = vector(0, m + q_degree);
  nrb_model->knv.m = m + q_degree;

  for(i = 0; i <= 3; i++)
    nrb_model->knu.U[i] = 0.0;
  for(i = nrb_model->net.n; i <= n+3; i++)
    nrb_model->knu.U[i] = 1.0;
  for(i = 4; i <= n-1; i++)
    nrb_model->knu.U[i] = (float)(i-3.0)/ (n-3.0);

  for(i = 0; i <= 3; i++)
    nrb_model->knv.U[i] = 0.0;
  for(i = nrb_model->net.m; i <= m+3; i++)
    nrb_model->knv.U[i] = 1.0;
  for(i = 4; i <= m-1; i++)
    nrb_model->knv.U[i] = (float)(i-3.0)/ (m-3.0);

  for(i = 0; i < m; i++)
    for(j = 0; j < n; j++)
      {
      nrb_model->net.Pw[j][i].x = 0.0;
      nrb_model->net.Pw[j][i].y = 0.0;
      nrb_model->net.Pw[j][i].z = 0.0;
      nrb_model->net.Pw[j][i].w = 0.0;
      }
}
/*-------------------------------------------------------------------------------*/   

/*----------------------------------------------------------------*/
void Find_miny(int m, SURFACE *nrb_surface, int *n)
/*----------------------------------------------------------------*/
/*----------------------------------------------------------------------
**  This subroutine is used to find the minimum value of y of a NURBS
**  surface at the specified m value.  The routine returns the n value
**  of the minimum
**----------------------------------------------------------------------
*/
{
  int j;
  float miny;
  
  miny = nrb_surface->net.Pw[0][m].y;
  *n = 0;

  for(j = 1; j < nrb_surface->net.n; j++)
    {
    if(nrb_surface->net.Pw[j][m].y < miny) 
      {
      miny = nrb_surface->net.Pw[j][m].y;
      *n = j;   
      }
    }
}

/*----------------------------------------------------------------*/
void Find_maxy(int m, SURFACE *nrb_surface, int *n)
/*----------------------------------------------------------------*/
/*----------------------------------------------------------------------
**  This subroutine is used to find the maximum value of y of a NURBS
**  surface at the specified m value.  The routine returns the n value
**  of the maximum
**----------------------------------------------------------------------
*/
{
  int j; 
  float maxy;
  
  maxy = nrb_surface->net.Pw[0][m].y;
  *n = 0;
    
  for(j = 1; j < nrb_surface->net.n; j++)
    {
    if(nrb_surface->net.Pw[j][m].y > maxy)
      {
      maxy = nrb_surface->net.Pw[j][m].y;
      *n = j;
      }
    }
}

/*----------------------------------------------------------------*/
void Find_minx(int m, SURFACE *nrb_surface, int *n)
/*----------------------------------------------------------------*/
/*----------------------------------------------------------------------
**  This subroutine is used to find the minimum value of x of a NURBS
**  surface at the specified m value.  The routine returns the n value
**  of the minimum
**----------------------------------------------------------------------
*/
{
  int j;
  float minx;
  
  minx = nrb_surface->net.Pw[0][m].x;
  *n = 0;
     
  for(j = 1; j < nrb_surface->net.n; j++) 
    {  
    if(nrb_surface->net.Pw[j][m].x < minx)
      {
      minx = nrb_surface->net.Pw[j][m].x;
      *n = j;
      }
    }
}


/*----------------------------------------------------------------*/
void Find_maxx(int m, SURFACE *nrb_surface, int *n)
/*----------------------------------------------------------------*/
/*----------------------------------------------------------------------
**  This subroutine is used to find the maximum value of x of a NURBS
**  surface at the specified m value.  The routine returns the n value
**  of the maximum
**----------------------------------------------------------------------
*/
{
  int j;
  float maxx;
  
  maxx = nrb_surface->net.Pw[0][m].x;
  *n = 0;
       
  for(j = 1; j < nrb_surface->net.n; j++) 
    {  
    if(nrb_surface->net.Pw[j][m].x > maxx)
      {
      maxx = nrb_surface->net.Pw[j][m].x;
      *n = j;
      }
    }
}

void CalculateCenter(SURFACE *nrb_surf, float *centerx, float *centery, float *centerz)
{
  int i, j;
       
  *centerx = 0.0; *centery = 0.0; *centerz = 0.0;
      
  for(i = 0; i < nrb_surf->net.m; i++)
    for(j = 0; j < nrb_surf->net.n; j++)
      {
      if(j != 1 && j != nrb_surf->net.n-1 && j != nrb_surf->net.n-2)
        {
        *centerx += nrb_surf->net.Pw[j][i].x;
        *centery += nrb_surf->net.Pw[j][i].y;
        *centerz += nrb_surf->net.Pw[j][i].z;
        }
      }

  *centerx /= ((nrb_surf->net.n-3.0) * nrb_surf->net.m);
  *centery /= ((nrb_surf->net.n-3.0) * nrb_surf->net.m);
  *centerz /= ((nrb_surf->net.n-3.0) * nrb_surf->net.m);
}

/*
void Read_torso_file(char *filename)
{
  FILE *fp;
  int i, j, k;
  int n, m;
  float temp, tx, ty, tz;
  char line[30];
  
  if((fp = fopen(filename, "r")) == NULL)
    Abort("Can not open torso nurbs datafile file");

  for(k = 0; k <= END_MODELS+12 ; k++)  { 
    fscanf(fp, "%s", line);

    fscanf(fp, "%i", &m); fscanf(fp, "%s", line);
    fscanf(fp, "%i", &n); fscanf(fp, "%s", line);

    initialize_structure(&nrb_model[k], n, m);

    fscanf(fp, "%s", line); fscanf(fp, "%s", line); fscanf(fp, "%s", line);
    for(i = 0; i <= n+p_degree; i++)
      {
      fscanf(fp, "%f", &temp);
      nrb_model[k].knu.U[i] = temp;
      }

    fscanf(fp, "%s", line); fscanf(fp, "%s", line); fscanf(fp, "%s", line);
    for(i = 0; i <= m+q_degree; i++)
      {
      fscanf(fp, "%f", &temp);
      nrb_model[k].knv.U[i] = temp;
      }
  
    fscanf(fp, "%s", line); fscanf(fp, "%s", line);
    for(i = 0; i < m; i++)
      for(j = 0; j< n; j++)  
        {
        fscanf(fp, "%f %f %f", &tx, &ty, &tz);

        tx += xoff;
        ty += yoff;
        tz += zoff;

        nrb_model[k].net.Pw[j][i].x = tx;
        nrb_model[k].net.Pw[j][i].y = ty;
        nrb_model[k].net.Pw[j][i].z = tz;
        nrb_model[k].net.Pw[j][i].w = 0.0;
        }
  }

  fscanf(fp, "%i", &gender_flag);
  fclose(fp);
}
*/

void Write(int intensity, int x, int y, float *out)
{
 if(x < 0)
   x = 0;
 if(x >= txdim)
   x = txdim -1;
 if(y < 0)
   y = 0;
 if(y >= tydim)
   y = tydim -1;
 
 if ( (0<=(y)) && ((y)<tydim) && (0<=(x)) && ((x)<txdim))
   out[(y)*txdim + (x)] = intensity;
}


void G_line(int x,int y,int x2,int y2, int intensity, float *out)
{
 int dx,dy,long_d,short_d;
 int d,add_dh,add_dl;
 register int inc_xh,inc_yh,inc_xl,inc_yl;
 register int i;

 /* Fix y value for screen */
 y = tydim - y;
 y2 = tydim - y2;

 dx=x2-x; dy=y2-y;                          /* ranges */

 if(dx<0){dx = -dx; inc_xh = -1; inc_xl = -1;}    /* making sure dx and dy >0 */
 else    {        inc_xh = 1;  inc_xl = 1; }    /* adjusting increments */
 if(dy<0){dy = -dy; inc_yh = -1; inc_yl = -1;}
 else    {        inc_yh = 1;  inc_yl = 1; }
          
 if(dx>dy){long_d=dx; short_d=dy; inc_yl=0;} /* long range,&making sure either */
 else     {long_d=dy; short_d=dx; inc_xl=0;} /* x or y is changed in L case */

 d=2*short_d-long_d;                        /* initial value of d */
 add_dl=2*short_d;                          /* d adjustment for H case */
 add_dh=2*short_d-2*long_d;                 /* d adjustment for L case */

 for(i=0;i<=long_d;i++)                     /* for all points in longer range */
 {
  Write(intensity,x,y, out);                        /* rendering */

  if(d>=0){x+=inc_xh; y+=inc_yh; d+=add_dh;}/* previous point was H type*/
  else    {x+=inc_xl; y+=inc_yl; d+=add_dl;}/* previous point was L type*/
 }
}


void lp_intersection(double *p0, double *p1, int z, int *result, int *flag)
{
  double t;

  *flag = 0;
  /*return with flag = 0 means no intersection*/
  /*return with flag = 1 means there is an intersection*/
  /*return with flag = 2 means line is on the plane*/

  result[0] = 0; result[1] = 0; result[2] = 0;

  if(p1[2]!= p0[2])
    {
    t = ( (double)z - p0[2]) / (p1[2] - p0[2]);
    if(t >= 0 && t <= 1)
      {
      *flag = 1;
      result[0] = floor( p0[0] + t * (p1[0] - p0[0]) +0.5);
      result[1] = floor( p0[1] + t * (p1[1] - p0[1]) +0.5);
      result[2] = z;    
      }
    }
  else if(p0[2] == z)
    *flag = 2;
}


void cross_product(double *u, double *v, double *result)
{
  result[0] = u[1]*v[2] - u[2]*v[1];
  result[1] = -(u[0]*v[2] - u[2]*v[0]);
  result[2] = u[0]*v[1] - u[1]*v[0];
}   

void Plane_eqn(double *p1, double *p2, double *p3, double *A, double *B, double *C, double *D)
{
  double v1[3], v2[3];
  double result[3];

  v1[0] = p2[0] - p1[0];
  v1[1] = p2[1] - p1[1];
  v1[2] = p2[2] - p1[2];

  v2[0] = p3[0] - p1[0];
  v2[1] = p3[1] - p1[1];
  v2[2] = p3[2] - p1[2];

  cross_product(v1, v2, result);

  *A = result[0];
  *B = result[1];
  *C = result[2];

  *D = -result[0] * p1[0] - result[1] * p1[1] - result[2] * p1[2];  
}


int Test_patch_z(double patch[4][4][3], int z)
{
  int i, j;
  double A, B, C, D;
  double denom;
  double t_error;
  double max_error;
  int flag = 0;

  Plane_eqn(patch[0][0], patch[0][3], patch[3][0], &A, &B, &C, &D); 
  denom = sqrt(A * A + B * B + C * C);
  if(denom == 0.0)
    {
    Plane_eqn(patch[0][3], patch[3][0], patch[3][3], &A, &B, &C, &D); 
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    {
    Plane_eqn(patch[0][0], patch[0][3], patch[3][3], &A, &B, &C, &D); 
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    {
    Plane_eqn(patch[0][0], patch[3][0], patch[3][3], &A, &B, &C, &D); 
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    {
    i = 0;
    while(denom == 0.0 && i <= 2)
      {
      j = 0;
      while(denom == 0.0 && j <= 2)
        {
        Plane_eqn(patch[i][j], patch[i][j+1], patch[i+1][j], &A, &B, &C, &D);
        j++;
        }
      i++;
      }  
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    {
    while(denom == 0.0 && i <= 2)
      {
      j = 0;
      while(denom == 0.0 && j <= 2)
        {
        Plane_eqn(patch[i][j], patch[i][j+1], patch[i+1][j+1], &A, &B, &C, &D);
        j++;
        }
      i++;
      }  
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    return 1;

  max_error = 0.0;
  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      {   
      if(denom != 0)
        t_error = fabs( (A * patch[i][j][0] + B * patch[i][j][1] + C * patch[i][j][2] + D) / denom);
      else
        t_error = 0.0;
      if(t_error > max_error)
        max_error = t_error;
      }

  if(max_error < 0.001 )
    return 1;
  else
    return 0;
  }

void hull_split_u(double P[4][4][3], double Q[4][4][3], double R[4][4][3])
{
  int i, iv;
  
  for(iv = 0; iv < 4; iv++)
    {
    for(i = 0; i < 3; i++)
      {
      Q[0][iv][i] = P[0][iv][i];
      Q[1][iv][i] = (P[0][iv][i]+P[1][iv][i]) / 2.0;
      Q[2][iv][i] = Q[1][iv][i] / 2.0 + (P[1][iv][i] + P[2][iv][i]) / 4.0;

      R[3][iv][i] = P[3][iv][i];
      R[2][iv][i] = (P[2][iv][i]+P[3][iv][i]) / 2.0;
      R[1][iv][i] = R[2][iv][i] / 2.0 + (P[1][iv][i] + P[2][iv][i]) / 4.0;

      Q[3][iv][i] = (Q[2][iv][i] + R[1][iv][i]) / 2.0;
      R[0][iv][i] = Q[3][iv][i];
      }
    }
}

void hull_split_v(double P[4][4][3], double Q[4][4][3], double R[4][4][3])
{
  int i, iu;
  
  for(iu = 0; iu < 4; iu++)
    {
    for(i = 0; i < 3; i++)
      {
      Q[iu][0][i] = P[iu][0][i];
      Q[iu][1][i] = (P[iu][0][i]+P[iu][1][i]) / 2.0;
      Q[iu][2][i] = Q[iu][1][i] / 2.0 + (P[iu][1][i] + P[iu][2][i]) / 4.0;

      R[iu][3][i] = P[iu][3][i];
      R[iu][2][i] = (P[iu][2][i]+P[iu][3][i]) / 2.0;
      R[iu][1][i] = R[iu][2][i] / 2.0 + (P[iu][1][i] + P[iu][2][i]) / 4.0;

      Q[iu][3][i] = (Q[iu][2][i] + R[iu][1][i]) / 2.0;
      R[iu][0][i] = Q[iu][3][i];
      }
    }
}

int count = 0;

void Render_triangle(TRIANGLE tri, int intensity, float *out, int z)
{
  int flag1, flag2, flag3, j;
  int int1[3], int2[3], int3[3];
  double v[3][3];
    
  for(j = 0; j < 3; j++)
    {
    v[j][0] = (tri.vertex[j].x*subvoxel_index) / (10 * pixel_width) + xoffset;
    v[j][1] = (tri.vertex[j].y*subvoxel_index) / (10 * pixel_width) + yoffset;
    v[j][2] = (tri.vertex[j].z*subvoxel_index) / (10 * slice_width) + zoffset;
    }
  
  /*Render triangle*/
  lp_intersection(v[0], v[2], z, int1, &flag1);
  lp_intersection(v[1], v[2], z, int2, &flag2);
  lp_intersection(v[0], v[1], z, int3, &flag3);
  
/*
  if(flag1 == 1 && flag2 == 1)
    G_line(int1[0], int1[1], int2[0], int2[1], intensity, out);
  else if(flag1 == 1 && flag3 == 1)
    G_line(int1[0], int1[1], int3[0], int3[1], intensity, out);
  else if(flag3 == 1 && flag2 == 1)
    G_line(int3[0], int3[1], int2[0], int2[1], intensity, out);
*/

  if(flag1 == 1 && flag2 == 1)
    G_line(int1[0], tydim-int1[1], int2[0], tydim-int2[1], intensity, out);
  else if(flag1 == 1 && flag3 == 1)
    G_line(int1[0], tydim-int1[1], int3[0], tydim-int3[1], intensity, out);
  else if(flag3 == 1 && flag2 == 1)
    G_line(int3[0], tydim-int3[1], int2[0], tydim-int2[1], intensity, out);
}

void Render_patch(double patch[4][4][3], int intensity, float *out, int z)
{
  int flag1, flag2, flag3;
  int int1[3], int2[3], int3[3];
                         
  /*Render 1st triangle*/
  lp_intersection(patch[3][0], patch[3][3], z, int1, &flag1);
  lp_intersection(patch[0][0], patch[3][3], z, int2, &flag2);
  lp_intersection(patch[3][0], patch[0][0], z, int3, &flag3);
  

  if(flag1 == 1 && flag2 == 1)
    G_line(int1[0], int1[1], int2[0], int2[1], intensity, out);
  else if(flag1 == 1 && flag3 == 1)
    G_line(int1[0], int1[1], int3[0], int3[1], intensity, out);
  else if(flag3 == 1 && flag2 == 1)
    G_line(int3[0], int3[1], int2[0], int2[1], intensity, out);

  /*Render 2nd triangle*/
  lp_intersection(patch[0][0], patch[3][3], z, int1, &flag1);
  lp_intersection(patch[0][3], patch[3][3], z, int2, &flag2);
  lp_intersection(patch[0][0], patch[0][3], z, int3, &flag3);

  if(flag1 == 1 && flag2 == 1)
    G_line(int1[0], int1[1], int2[0], int2[1], intensity, out);
  else if(flag1 == 1 && flag3 == 1)
    G_line(int1[0], int1[1], int3[0], int3[1], intensity, out);
  else if(flag3 == 1 && flag2 == 1)
    G_line(int3[0], int3[1], int2[0], int2[1], intensity, out); 
}


void Subdivide_patch(double patch[4][4][3], double ul_patch[4][4][3], double ur_patch[4][4][3], double dl_patch[4][4][3],
double dr_patch[4][4][3])
{
  double up_patch[4][4][3], down_patch[4][4][3];
  int i, j, k;

  hull_split_u(patch, up_patch, down_patch); 

  hull_split_v(up_patch, ul_patch, ur_patch); 
  hull_split_v(down_patch, dl_patch, dr_patch); 
}


int Test_extents_z(double patch[4][4][3], int z)
{
  int i, j;
  double minz, maxz;
  double z_test;

  z_test = z;

  minz = patch[0][0][2];
  maxz = minz;

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      {
       if(patch[i][j][2] < minz)
         minz = patch[i][j][2];
       if(patch[i][j][2] > maxz)
         maxz = patch[i][j][2];
      }

  if(z_test >= minz && z_test <= maxz)
    return 1;
  else
   return 0;
}

void Intersect_bez_z(double patch[4][4][3], float *out, int intensity, int z)
{
  double ul_patch[4][4][3], ur_patch[4][4][3], dl_patch[4][4][3], dr_patch[4][4][3];

  if( Test_patch_z(patch, z) )
    Render_patch(patch, intensity, out, z);
  else
    {
    Subdivide_patch(patch, ul_patch, ur_patch, dl_patch, dr_patch);
    if(Test_extents_z(ul_patch, z))
      Intersect_bez_z(ul_patch, out, intensity, z);
    if(Test_extents_z(ur_patch, z))
      Intersect_bez_z(ur_patch, out, intensity, z);  
    if(Test_extents_z(dl_patch, z))
      Intersect_bez_z(dl_patch, out, intensity, z);  
    if(Test_extents_z(dr_patch, z))
      Intersect_bez_z(dr_patch, out, intensity, z);  
    }
}


void Render_Bezier(BEZIER_MODEL *bez_model, float *out, int intensity, int z)
{
  int i, j, k, l;
  float temp1, temp2, temp3, temp4;
  double patch[4][4][3];

  for(i = 0; i < bez_model->num_patches; i++)
    {
    if(z >= bez_model->patches[i].slice_minz && z <= bez_model->patches[i].slice_maxz)
      Intersect_bez_z(bez_model->patches[i].slice_points, out, intensity, z);
    }
}

int Check_Y_Boundary(int x, int y, int threshold, float *out)
{
  int bound1=0, bound2=0;
  int i, j;
  unsigned long int index;

  i = y;
  while(!bound1 && i < tydim)
    {
    index = x + i * txdim;
    if( (int)out[index] == threshold)
      bound1 = 1;
    i++;
    }
 
  i = y;
  while(!bound2 && i >= 0)
    {
    index = x+ i * txdim;     
    if( (int)out[index] == threshold)
      bound2 = 1;
    i--;
    }

  return(bound1&&bound2);
}

void Dabs(double num, double *result)
{
  if (num > 0.0)
    *result = num;
  else
   *result = -num;
}


/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/
template<typename FloatType, typename FloatArrayType>
void Rotate_Y(FloatType cosRot, FloatType sinRot, FloatArrayType p)
{
	// Rotate about y
	const FloatType x = p[0] * cosRot - p[2] * sinRot;
	const FloatType z = p[0] * sinRot + p[2] * cosRot;

	p[2] = z;
	p[0] = x;
}
template<typename FloatType, typename FloatArrayType>
void Rotate_Y(FloatType y_rot, FloatArrayType p)
/*---------------------------------------------------------------------------------------
**  This subroutine is used to rotate a point (p) about the y axis.
**  y_rot = angle to rotate about the y-axis
**---------------------------------------------------------------------------------------
*/
{
  const FloatType angle_rad = y_rot * (PI / 180.0f);
  const FloatType cosine = cos(angle_rad);
  const FloatType sine = sin(angle_rad);

  // Rotate about y
  Rotate_Y(cosine, sine, p);
}
/*---------------------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------*/
template<typename FloatType, typename FloatArrayType>
void Rotate_Z(FloatType cosRot, FloatType sinRot, FloatArrayType p)
/*---------------------------------------------------------------------------------------
**  This subroutine is used to rotate a point (p) about the z axis.
**  z_rot = angle to rotate about the z-axis
**---------------------------------------------------------------------------------------
*/
{
	// Rotate about z
	const FloatType x = p[0] * cosRot + p[1] * sinRot;
	const FloatType y = -p[0] * sinRot + p[1] * cosRot;

	p[1] = y;
	p[0] = x;
}
template<typename FloatType, typename FloatArrayType>
void Rotate_Z(FloatType z_rot, FloatArrayType p)
/*---------------------------------------------------------------------------------------
**  This subroutine is used to rotate a point (p) about the z axis.
**  z_rot = angle to rotate about the z-axis
**---------------------------------------------------------------------------------------
*/
{
  const FloatType angle_rad = z_rot * (PI / 180.0f);
  const FloatType cosine = cos(angle_rad);
  const FloatType sine = sin(angle_rad);
  
  // Rotate about z
  Rotate_Z(cosine, sine, p);
}
/*---------------------------------------------------------------------------------------*/


void Fill(int intensity, float *out)
{   
   int x, y;
   unsigned long index;
   int diff;
   int *edge;
   int edge_counter;
   int counter = 0;
   int flag = 0;
   int i;

   edge = ivector(0, XDIM*subvoxel_index-1);
   for(i = 0; i < XDIM*subvoxel_index; i++)
     edge[i] = -1;

   for(y = 0; y < tydim; y++) {  
     edge_counter = 0;
     for(x = 0; x < txdim; x++) {
       index = x+ y * txdim;     
       diff = (int)out[index];
       if( diff == intensity) 
         {
         edge[edge_counter] = x;
         edge_counter++;
         if(edge_counter == XDIM*subvoxel_index)
           {
           printf("\nERROR: Edge_counter = %i", edge_counter);
           exit(1);
           }
         }
     }

     while(edge_counter!=0)
       {
       flag = 0;
       if(edge_counter-2 >=0)
         {
         if(edge[edge_counter-1] != -1 && edge[edge_counter-2] != -1)
           {
           if(edge[edge_counter-1] - edge[edge_counter-2] > 1)
             flag = Check_Y_Boundary( (edge[edge_counter-2]+edge[edge_counter-1])/2, y,intensity,out);

           if(flag)
             {           
             G_line(edge[edge_counter-2], tydim - y, edge[edge_counter-1], tydim - y, intensity, out);
             flag = 0;
             }
           }
         }
       edge_counter--;
       }
   }

  free_ivector(edge, 0, XDIM*subvoxel_index-1);
}



/*-------------------------------------------------------------------------*/
int Test_patch(double patch[4][4][3], double tol)
/*-------------------------------------------------------------------------*/
/* This function tests a Bezier patch to see if it is flat and can be      */
/* approximated by a rectangle                                             */
/*-------------------------------------------------------------------------*/
{
  int i, j;
  double A, B, C, D;
  double denom;
  double t_error;
  double max_error;

  Plane_eqn(patch[0][0], patch[0][3], patch[3][0], &A, &B, &C, &D); 
 
  denom = sqrt(A * A + B * B + C * C);
  if(denom == 0.0)
    {
    Plane_eqn(patch[0][3], patch[3][0], patch[3][3], &A, &B, &C, &D); 
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    {
    Plane_eqn(patch[0][0], patch[0][3], patch[3][3], &A, &B, &C, &D); 
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    {
    Plane_eqn(patch[0][0], patch[3][0], patch[3][3], &A, &B, &C, &D); 
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    {
    i = 0;
    while(denom == 0.0 && i <= 2)
      {
      j = 0;
      while(denom == 0.0 && j <= 2)
        {
        Plane_eqn(patch[i][j], patch[i][j+1], patch[i+1][j], &A, &B, &C, &D);
        j++;
        }
      i++;
      }  
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    {
    while(denom == 0.0 && i <= 2)
      {
      j = 0;
      while(denom == 0.0 && j <= 2)
        {
        Plane_eqn(patch[i][j], patch[i][j+1], patch[i+1][j+1], &A, &B, &C, &D);
        j++;
        }
      i++;
      }  
    denom = sqrt(A * A + B * B + C * C);
    }
  if(denom == 0.0)
    return 1;

  max_error = -1.0;

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      {   
      if(denom != 0.0)
        t_error = fabs( (A * patch[i][j][0] + B * patch[i][j][1] + C * patch[i][j][2] + D) / denom);
      else
        t_error = 0.0;
      if(t_error > max_error)
        max_error = t_error;
      }

  if(max_error < tol)
    return 1;
  else
    return 0;
  }

int Test_extents(float xl, float xh, float yl, float yh, float zl, float zh, float line_origin[3], float line_vector[3])
{
  int flag = 0;
  int i, j;

  float x0,y0,z0;
  float xd,yd,zd;

  float t1, t2, tnear = -100000, tfar = 100000;
  float temp;

  x0 = line_origin[0]; y0 = line_origin[1]; z0 = line_origin[2];
  xd = line_vector[0]; yd = line_vector[1]; zd = line_vector[2];

  /* X PLANES */
  if(xd == 0)
    {
    if(x0 < xl || x0 > xh)
      return 0;
    }
  else
    {
    t1 = (xl-x0)/xd;
    t2 = (xh-x0)/xd;

    if(t1 > t2)
      {
      temp = t1;
      t1 = t2;
      t2 = temp;
      }
    if(t1 > tnear)
      tnear = t1;
    if(t2 < tfar)
      tfar = t2;

    flag = 1;
    if(tnear > tfar)
      return 0;
    if(tfar < 0)
      return 0;
    }

  /* Y PLANES */
  if(yd == 0)
    {
    if(y0 < yl || y0 > yh)
      return 0;
    }
  else
    {
    t1 = (yl-y0)/yd;
    t2 = (yh-y0)/yd;

    if(t1 > t2)
      {
      temp = t1;
      t1 = t2;
      t2 = temp;
      }
    if(t1 > tnear)
      tnear = t1;
    if(t2 < tfar)
      tfar = t2;

    flag = 1;
    if(tnear > tfar)
      return 0;
    if(tfar < 0)
      return 0;
    }

  /* Z PLANES */
  if(zd == 0)
    {
    if(z0 < zl || z0 > zh)
      return 0;
    }
  else
    {
    t1 = (zl-z0)/zd;
    t2 = (zh-z0)/zd;

    if(t1 > t2)
      {
      temp = t1;
      t1 = t2;
      t2 = temp;
      }
    if(t1 > tnear)
      tnear = t1;
    if(t2 < tfar)
      tfar = t2;

    flag = 1;
    if(tnear > tfar)
      return 0;

    if(tfar < 0)
      return 0;
    }

  return 1;
}

int Test_extents2(double patch[4][4][3], float line_origin[3], float line_vector[3])
{
  int flag = 0;
  float xl,yl,zl;
  float xh,yh,zh;
  int i, j;

  float x0,y0,z0;
  float xd,yd,zd;

  float t1, t2, tnear = -100000, tfar = 100000;
  float temp;

  zl = patch[0][0][2];
  zh = zl;

  yl = patch[0][0][1];
  yh = yl;

  xl = patch[0][0][0];
  xh = xl;

  for(i = 0; i < 4; i++)
    for(j = 0; j < 4; j++)
      {
       if(patch[i][j][2] < zl)
         zl = patch[i][j][2];
       if(patch[i][j][2] > zh)
         zh = patch[i][j][2];

       if(patch[i][j][1] < yl)
         yl = patch[i][j][1];
       if(patch[i][j][1] > yh)
         yh = patch[i][j][1];

       if(patch[i][j][0] < xl)
         xl = patch[i][j][0];
       if(patch[i][j][0] > xh)
         xh = patch[i][j][0];
      }

  x0 = line_origin[0]; y0 = line_origin[1]; z0 = line_origin[2];
  xd = line_vector[0]; yd = line_vector[1]; zd = line_vector[2];

  /* X PLANES */
  if(xd == 0)
    {
    if(x0 < xl || x0 > xh)
      return 0;
    }
  else
    {
    t1 = (xl-x0)/xd;
    t2 = (xh-x0)/xd;

    if(t1 > t2)
      {
      temp = t1;
      t1 = t2;
      t2 = temp;
      }
    if(t1 > tnear)
      tnear = t1;
    if(t2 < tfar)
      tfar = t2;

    flag = 1;
    if(tnear > tfar)
      return 0;
    if(tfar < 0)
      return 0;
    }

  /* Y PLANES */
  if(yd == 0)
    {
    if(y0 < yl || y0 > yh)
      return 0;
    }
  else
    {
    t1 = (yl-y0)/yd;
    t2 = (yh-y0)/yd;

    if(t1 > t2)
      {
      temp = t1;
      t1 = t2;
      t2 = temp;
      }
    if(t1 > tnear)
      tnear = t1;
    if(t2 < tfar)
      tfar = t2;

    flag = 1;
    if(tnear > tfar)
      return 0;
    if(tfar < 0)
      return 0;
    }

  /* Z PLANES */
  if(zd == 0)
    {
    if(z0 < zl || z0 > zh)
      return 0;
    }
  else
    {
    t1 = (zl-z0)/zd;
    t2 = (zh-z0)/zd;

    if(t1 > t2)
      {
      temp = t1;
      t1 = t2;
      t2 = temp;
      }
    if(t1 > tnear)
      tnear = t1;
    if(t2 < tfar)
      tfar = t2;

    flag = 1;
    if(tnear > tfar)
      return 0;

    if(tfar < 0)
      return 0;
    }

  return 1;
}

void Check_difference(double xint, double *min_error, XP_ARRAY *int_points)
{
  int i;
  double tmp_error;
  
  *min_error = 1.0;

  for(i = 0; i < int_points->length; i++)
    {
    tmp_error = fabs(xint - int_points->xp[i].x);
    if(tmp_error < *min_error)
      *min_error = tmp_error;
    }
}


void Intersect_bez(double patch[4][4][3], int organ_id, XP_ARRAY *int_points, double tol, float line_origin[3],
float line_vector[3])
/* ------------------------------------------------------------------------------------ */
/* This routine iteratively breaks down surface patches until it finds one that is flat */
/* and is intersected by the projection ray.  Once it finds a flat patch, it calculates */
/* the intersection point                                                               */
/* ------------------------------------------------------------------------------------ */
{
  double ul_patch[4][4][3], ur_patch[4][4][3], dl_patch[4][4][3], dr_patch[4][4][3];
  double A,B,C,D;
  double numerator, denominator;
  double xint, min_error;
  int count, count2;
  int i, j;
  int flag;

  if(Test_patch(patch, tol))
    {
    Plane_eqn(patch[0][0], patch[0][3], patch[3][0], &A, &B, &C, &D);

    numerator   = A*line_origin[0] + B*line_origin[1] + C*line_origin[2] + D;
    denominator = -A*line_vector[0] - B*line_vector[1] - C*line_vector[2];

    if(denominator != 0)
      xint = numerator/denominator;
    else
      xint = 0;

    if(int_points->length == 0 && xint > 0)
      {
      int_points->length = 1;
      int_points->xp[0].x = xint;
      int_points->xp[0].organ_id = organ_id;
      }
    else if(xint > 0)
      {
      Check_difference(xint, &min_error, int_points);
      if(min_error >= 1.0)
        {
        count = 0;
        flag = 1;
        while(count < int_points->length && flag)
          if(xint > int_points->xp[count].x)
            count++;
          else
            flag = 0;
        for(count2 = int_points->length; count2 > count; count2--)
          {
          int_points->xp[count2].x = int_points->xp[count2-1].x;
          int_points->xp[count2].organ_id = int_points->xp[count2-1].organ_id;
          }
        int_points->xp[count].x = xint;
        int_points->xp[count].organ_id = organ_id;
        int_points->length++;
        }
      }
    }
  else
    {
    Subdivide_patch(patch, ul_patch, ur_patch, dl_patch, dr_patch);
    if(Test_extents2(ul_patch, line_origin, line_vector))
      Intersect_bez(ul_patch, organ_id, int_points, tol, line_origin, line_vector);
    if(Test_extents2(ur_patch, line_origin, line_vector))
      Intersect_bez(ur_patch, organ_id, int_points, tol, line_origin, line_vector);
    if(Test_extents2(dl_patch, line_origin, line_vector))
      Intersect_bez(dl_patch, organ_id, int_points, tol, line_origin, line_vector);
    if(Test_extents2(dr_patch, line_origin, line_vector))
      Intersect_bez(dr_patch, organ_id, int_points, tol, line_origin, line_vector);
    }
}

void Find_Intersections(BEZIER_MODEL *bez_model, int organ_id, float line_origin[3], float line_vector[3], XP_ARRAY
*int_points, double tol)
{
  int i;

  for(i = 0; i < bez_model->num_patches; i++)
    {
    if(Test_extents(bez_model->patches[i].minx, bez_model->patches[i].maxx,
                    bez_model->patches[i].miny, bez_model->patches[i].maxy,
                    bez_model->patches[i].minz, bez_model->patches[i].maxz,
                    line_origin, line_vector))
     Intersect_bez(bez_model->patches[i].cntrl_points, organ_id, int_points, tol, line_origin, line_vector);
    }
}

float magnitude(float *u)
{
  return sqrt(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]);
}
  
void normalize(float *u)
{
  float mag = magnitude(u);   
    
  u[1] /= mag;
  u[2] /= mag;
  u[3] /= mag;
}

float dot_product(float *u, float *v)
{
  float temp;
  temp = (u[1]*v[1] + u[2]*v[2] + u[3]*v[3]);
  if(temp <= -1.0)
    temp = -1.0;
  if(temp >= 1.0)
    temp = 1.0;
  return temp;
}
   
void GetVolumes(float *out_pixel)
{
  int ID;
  float vol[END_MODELS+1];

  for(unsigned index = 0; index <= END_MODELS; index++)
    vol[index] = 0.0;

  for(unsigned index = 0; index < XDIM*YDIM*ZDIM_OUTPUT; index++)
  {
    ID = out_pixel[index];
    vol[ID]++;
  }

printf("\nVolume of superior frontal gyrus right                                    %f", vol[10]);
printf("\nVolume of superior frontal gyrus left                                     %f", vol[70]);
printf("\nVolume of middle frontal gyrus right                                      %f", vol[2]);
printf("\nVolume of middle frontal gyrus left                                       %f", vol[50]);
printf("\nVolume of inferior frontal gyrus right                                    %f", vol[75]);
printf("\nVolume of inferior frontal gyrus left                                     %f", vol[15]);
printf("\nVolume of precentral gyrus right                                          %f", vol[5]);
printf("\nVolume of precentral gyrus left                                           %f", vol[80]);
printf("\nVolume of lateral front-orbital gyrus right                               %f", vol[6]);
printf("\nVolume of lateral front-orbital gyrus left                                %f", vol[90]);
printf("\nVolume of medial front-orbital gyrus right                                %f", vol[1]);
printf("\nVolume of medial front-orbital gyrus left                                 %f", vol[85]);
printf("\nVolume of cingulate region right                                          %f", vol[7]);
printf("\nVolume of cingulate region left                                           %f", vol[27]);
printf("\nVolume of medial frontal gyrus right                                      %f", vol[114]);
printf("\nVolume of medial frontal gyrus left                                       %f", vol[9]);
printf("\nVolume of superior parietal lobule right                                  %f", vol[88]);
printf("\nVolume of superior parietal lobule left                                   %f", vol[52]);
printf("\nVolume of supramarginal gyrus right                                       %f", vol[60]);
printf("\nVolume of supramarginal gyrus left                                        %f", vol[41]);
printf("\nVolume of angular gyrus right                                             %f", vol[19]);
printf("\nVolume of angular gyrus left                                              %f", vol[159]);
printf("\nVolume of precuneus right                                                 %f", vol[32]);
printf("\nVolume of precuneus left                                                  %f", vol[56]);
printf("\nVolume of postcentral gyrus right                                         %f", vol[110]);
printf("\nVolume of postcentral gyrus left                                          %f", vol[74]);
printf("\nVolume of superior temporal gyrus right                                   %f", vol[145]);
printf("\nVolume of superior temporal gyrus left                                    %f", vol[61]);
printf("\nVolume of middle temporal gyrus right                                     %f", vol[130]);
printf("\nVolume of middle temporal gyrus left                                      %f", vol[64]);
printf("\nVolume of inferior temporal gyrus right                                   %f", vol[140]);
printf("\nVolume of inferior temporal gyrus left                                    %f", vol[164]);
printf("\nVolume of uncus right                                                     %f", vol[26]);
printf("\nVolume of uncus left                                                      %f", vol[62]);
printf("\nVolume of medial occipitotemporal gyrus right                             %f", vol[165]);
printf("\nVolume of medial occipitotemporal gyrus left                              %f", vol[119]);
printf("\nVolume of lateral occipitotemporal gyrus right                            %f", vol[99]);
printf("\nVolume of lateral occipitotemporal gyrus left                             %f", vol[196]);
printf("\nVolume of amygdala right                                                  %f", vol[139]);
printf("\nVolume of amygdala left                                                   %f", vol[118]);
printf("\nVolume of parahippocampal gyrus right                                     %f", vol[125]);
printf("\nVolume of parahippocampal gyrus left                                      %f", vol[18]);
printf("\nVolume of occipital pole right                                            %f", vol[132]);
printf("\nVolume of occipital pole left                                             %f", vol[251]);
printf("\nVolume of superior occipital gyrus right                                  %f", vol[38]);
printf("\nVolume of superior occipital gyrus left                                   %f", vol[98]);
printf("\nVolume of middle occipital gyrus right                                    %f", vol[63]);
printf("\nVolume of middle occipital gyrus left                                     %f", vol[154]);
printf("\nVolume of inferior occipital gyrus right                                  %f", vol[97]);
printf("\nVolume of inferior occipital gyrus left                                   %f", vol[37]);
printf("\nVolume of cuneus right                                                    %f", vol[175]);
printf("\nVolume of cuneus left                                                     %f", vol[54]);
printf("\nVolume of lingual gyrus right                                             %f", vol[112]);
printf("\nVolume of lingual gyrus left                                              %f", vol[69]);
printf("\nVolume of insula right                                                    %f", vol[4]);
printf("\nVolume of insula left                                                     %f", vol[108]);
printf("\nVolume of caudate nucleus right                                           %f", vol[53]);
printf("\nVolume of caudate nucleus left                                            %f", vol[39]);
printf("\nVolume of putamen right                                                   %f", vol[16]);
printf("\nVolume of putamen left                                                    %f", vol[14]);
printf("\nVolume of globus palladus right                                           %f", vol[11]);
printf("\nVolume of globus palladus left                                            %f", vol[12]);
printf("\nVolume of thalamus right                                                  %f", vol[203]);
printf("\nVolume of thalamus left                                                   %f", vol[102]);
printf("\nVolume of corpus callosum                                                 %f", vol[133]);
printf("\nVolume of subarachnoid cerebro-spinal fluid                               %f", vol[255]);
printf("\nVolume of third ventricle                                                 %f", vol[232]);
printf("\nVolume of fourth ventricle                                                %f", vol[233]);
printf("\nVolume of lateral ventricle right                                         %f", vol[8]);
printf("\nVolume of lateral ventricle left                                          %f", vol[3]);
printf("\nVolume of brain stem                                                      %f", vol[20]);
printf("\nVolume of frontal lobe WM right                                           %f", vol[17]);
printf("\nVolume of frontal lobe WM left                                            %f", vol[30]);
printf("\nVolume of occipital lobe WM right                                         %f", vol[45]);
printf("\nVolume of occipital lobe WM left                                          %f", vol[73]);
printf("\nVolume of parietal lobe WM right                                          %f", vol[105]);
printf("\nVolume of parietal lobe WM left                                           %f", vol[57]);
printf("\nVolume of temporal lobe WM right                                          %f", vol[59]);
printf("\nVolume of temporal lobe WM left                                           %f", vol[83]);
printf("\nVolume of scalp                                                           %f", vol[21]);
printf("\nVolume of skull                                                           %f", vol[28]);
printf("\nVolume of hippocampal formation right                                     %f", vol[36]);
printf("\nVolume of hippocampal formation left                                      %f", vol[101]);
printf("\nVolume of nucleus accumbens right                                         %f", vol[25]);
printf("\nVolume of nucleus accumbens left                                          %f", vol[72]);
printf("\nVolume of fornix right                                                    %f", vol[254]);
printf("\nVolume of fornix left                                                     %f", vol[29]);
printf("\nVolume of posterior limb of internal capsule inc. cerebral peduncle right %f", vol[35]);
printf("\nVolume of posterior limb of internal capsule inc. cerebral peduncle left  %f", vol[34]);
printf("\nVolume of subthalamic nucleus right                                       %f", vol[23]);
printf("\nVolume of subthalamic nucleus left                                        %f", vol[33]);
printf("\nVolume of anterior limb of internal capsule right                         %f", vol[128]);
printf("\nVolume of anterior limb of internal capsule left                          %f", vol[43]);
printf("\nVolume of cerebellum right                                                %f", vol[76]);
printf("\nVolume of cerebellum left                                                 %f", vol[67]);
}


/*------------------------Triangle Routines-----------------------*/
void Calc_extents_tri(TRI_MODEL *tmodel)
{
  int i, j;

  tmodel->min_x = 10000;
  tmodel->min_y = 10000;
  tmodel->min_z = 10000;

  tmodel->max_x = -10000;
  tmodel->max_y = -10000;
  tmodel->max_z = -10000;

  for(i = 0; i < tmodel->num_tris; i++)
    for(j = 0; j < 3; j++)
      {
      if(tmodel->tris[i].vertex[j].x < tmodel->min_x)
        tmodel->min_x = tmodel->tris[i].vertex[j].x;
      if(tmodel->tris[i].vertex[j].y < tmodel->min_y)
        tmodel->min_y = tmodel->tris[i].vertex[j].y;
      if(tmodel->tris[i].vertex[j].z < tmodel->min_z)
        tmodel->min_z = tmodel->tris[i].vertex[j].z;

      if(tmodel->tris[i].vertex[j].x > tmodel->max_x)
        tmodel->max_x = tmodel->tris[i].vertex[j].x;
      if(tmodel->tris[i].vertex[j].y > tmodel->max_y)
        tmodel->max_y = tmodel->tris[i].vertex[j].y;
      if(tmodel->tris[i].vertex[j].z > tmodel->max_z)
        tmodel->max_z = tmodel->tris[i].vertex[j].z;
      }

 tmodel->max_z = floor( (tmodel->max_z) / (slice_width/subvoxel_index * 10.0) + 0.5);
 tmodel->min_z = floor( (tmodel->min_z) / (slice_width/subvoxel_index * 10.0) + 0.5);
}

void MoveTriangleModel(TRI_MODEL *tmodel, float x_trans, float y_trans, float z_trans, float x_rot, float y_rot, float z_rot, float *tx, float *ty, float *tz)
{
  float p[4];
  int i, j;
  float R[3][3];
  float X, Y, Z;
  float newx, newy, newz;

  X = x_rot*PI/180.0;
  Y = y_rot*PI/180.0;
  Z = z_rot*PI/180.0;

  R[0][0] = cos(Y)*cos(Z); R[0][1] = sin(X)*sin(Y)*cos(Z) - cos(X)*sin(Z); R[0][2] = cos(X)*sin(Y)*cos(Z) + sin(X)*sin(Z);
  R[1][0] = cos(Y)*sin(Z); R[1][1] = sin(X)*sin(Y)*sin(Z) + cos(X)*cos(Z); R[1][2] = cos(X)*sin(Y)*sin(Z) - sin(X)*cos(Z);
  R[2][0] = -sin(Y);       R[2][1] = cos(Y)*sin(X);                        R[2][2] = cos(Y)*cos(X);


  for(i = 0; i < tmodel->num_tris; i++)
    {
    for(j = 0; j < 3; j++)
      {
      p[1] = tmodel->tris[i].vertex[j].x;
      p[2] = tmodel->tris[i].vertex[j].y;
      p[3] = tmodel->tris[i].vertex[j].z;

      p[1] -= *tx;
      p[2] -= *ty;
      p[3] -= *tz;

      newx = R[0][0]*p[1] + R[0][1]*p[2] + R[0][2]*p[3] + x_trans + *tx;
      newy = R[1][0]*p[1] + R[1][1]*p[2] + R[1][2]*p[3] + y_trans + *ty;
      newz = R[2][0]*p[1] + R[2][1]*p[2] + R[2][2]*p[3] + z_trans + *tz;

      tmodel->tris[i].vertex[j].x = newx;
      tmodel->tris[i].vertex[j].y = newy;
      tmodel->tris[i].vertex[j].z = newz;
      }
    }
  
  *tx += x_trans; *ty += y_trans; *tz += z_trans;

  Calc_extents_tri(tmodel);
}

void Read_Triangle_Model(char *filename, TRI_MODEL *tmodel)
{
  FILE *fp;
  int i, j, num;
  float x, y, z;

  /* Open input file */
  if( (fp = fopen(filename, "r")) == NULL )
    {
    tmodel->num_tris = 0;
    return;
    }

  fscanf(fp, "%i", &num);
  tmodel->num_tris = num;
  tmodel->tris = tri_vector(0, num-1);

  for(i = 0; i < num; i++)
    {
    for(j = 0; j < 3; j++)
       {
       fscanf(fp, "%f %f %f", &x, &y, &z);
       tmodel->tris[i].vertex[j].x = x + xoff;
       tmodel->tris[i].vertex[j].y = y + yoff;
       tmodel->tris[i].vertex[j].z = z + zoff;
       }
    }

  Calc_extents_tri(tmodel);
  fclose(fp);
}

void Output_Triangle_Model(char *filename, TRI_MODEL *tmodel)
{
  FILE *fp;
  int i, j, num;
  float x, y, z;

  /* Open input file */
  if( (fp = fopen(filename, "w")) == NULL )
    return;

  fprintf(fp, "Object1");

  for(i = 0; i < tmodel->num_tris; i++)
    {
    fprintf(fp, "\n");
    for(j = 0; j < 3; j++)
       fprintf(fp, "%f %f %f ", tmodel->tris[i].vertex[j].x, tmodel->tris[i].vertex[j].y, tmodel->tris[i].vertex[j].z);
    }

  fclose(fp);
}


void Read_Triangle_ModelVTK(char *filename, TRI_MODEL *tmodel)
{
  FILE *fp;
  int i, j, tri_count, num_tris, num_points, num_strips, num;
  float x, y, z;
  int fields[5000];
  POINT *P;
  char line[50];

  /* Open input file */
  if( (fp = fopen(filename, "r")) == NULL )
    {
    tmodel->num_tris = 0;
    return;
    }
  
  NEXTLINE(fp);
  NEXTLINE(fp);
  NEXTLINE(fp);
  NEXTLINE(fp);
  fscanf(fp, "%s %i %s", line, &num_points, line);
  P = p_vector(0, num_points-1);
  for(i = 0; i < num_points; i++)
    fscanf(fp, "%f %f %f", &P[i].x, &P[i].y, &P[i].z);

  fscanf(fp, "%s %i %s", line, &num_strips, line);
  num_tris = 0;
  for(i = 0; i < num_strips; i++)
    {
    fscanf(fp, "%i", &num);
    if(num >= 5000)
      {
      printf("\nNeed to increase size for triangle strips");
      exit(1);
      }
    for(j = 0; j < num; j++)
      fscanf(fp, "%i", &fields[j]);

    num_tris += num-2;
    }
  fclose(fp);

  tmodel->num_tris = num_tris;
  tmodel->tris = tri_vector(0, num_tris-1);

  /* Open input file */
  if( (fp = fopen(filename, "r")) == NULL )
    {
    tmodel->num_tris = 0;
    return;
    }

  NEXTLINE(fp);
  NEXTLINE(fp);
  NEXTLINE(fp);
  NEXTLINE(fp);
  fscanf(fp, "%s %i %s", line, &num_points, line);
  for(i = 0; i < num_points; i++)
    fscanf(fp, "%f %f %f", &P[i].x, &P[i].y, &P[i].z);

  fscanf(fp, "%s %i %s", line, &num_strips, line);
  tri_count = 0;
  for(i = 0; i < num_strips; i++)
    {
    fscanf(fp, "%i", &num);
    if(num >= 5000)
      {
      printf("\nNeed to increase size for triangle strips");
      exit(1);
      }
    for(j = 0; j < num; j++)
      fscanf(fp, "%i", &fields[j]);

    for(j = 0; j < num-2; j++)
      {
      tmodel->tris[tri_count].vertex[0].x = P[fields[j]].x + xoff;
      tmodel->tris[tri_count].vertex[0].y = P[fields[j]].y + yoff;
      tmodel->tris[tri_count].vertex[0].z = P[fields[j]].z + zoff;

      tmodel->tris[tri_count].vertex[1].x = P[fields[j+1]].x + xoff;
      tmodel->tris[tri_count].vertex[1].y = P[fields[j+1]].y + yoff;
      tmodel->tris[tri_count].vertex[1].z = P[fields[j+1]].z + zoff;

      tmodel->tris[tri_count].vertex[2].x = P[fields[j+2]].x + xoff;
      tmodel->tris[tri_count].vertex[2].y = P[fields[j+2]].y + yoff;
      tmodel->tris[tri_count].vertex[2].z = P[fields[j+2]].z + zoff;
      tri_count++;
      }
    }

  free_pvector(P, 0, num_points-1);
  Calc_extents_tri(tmodel);
  fclose(fp);
}


void Translate_tri(TRI_MODEL *tmodel, float tx, float ty, float tz)
{
  int i, j;

  for(i = 0; i < tmodel->num_tris; i++)
    for(j = 0; j < 3; j++)
      {
      tmodel->tris[i].vertex[j].x += tx;
      tmodel->tris[i].vertex[j].y += ty;
      tmodel->tris[i].vertex[j].z += tz;
      }
  Calc_extents_tri(tmodel);
}

int Test_extents_tri_z(TRIANGLE tri, int z)
{
  int j;
  double minz, maxz;

  minz = tri.vertex[0].z;
  maxz = minz;

  for(j = 1; j < 3; j++)
    {
     if(tri.vertex[j].z < minz)
       minz = tri.vertex[j].z;
     if(tri.vertex[j].z > maxz)
       maxz = tri.vertex[j].z;
    }

  minz = (minz*subvoxel_index) / (10 * slice_width) + zoffset;
  maxz = (maxz*subvoxel_index) / (10 * slice_width) + zoffset;

  if(z >= minz-1 && z <= maxz+1)
    return 1;
  else
   return 0;
}

int Test_extents_tri(TRIANGLE T)
/* Function tests the extents of a Triangle */
{
  int i, j;
  double minz, maxz;
  double miny, maxy;
  double minx, maxx;

  minz = T.vertex[0].z;
  maxz = minz;

  miny = T.vertex[0].y;
  maxy = miny;

  minx = T.vertex[0].x;
  maxx = minx;

  for(i = 1; i < 3; i++)
    {
    if(T.vertex[i].z < minz)
      minz = T.vertex[i].z;
    if(T.vertex[i].z > maxz)
      maxz = T.vertex[i].z;

    if(T.vertex[i].y < miny)
      miny = T.vertex[i].y;
    if(T.vertex[i].y > maxy)
      maxy = T.vertex[i].y;

    if(T.vertex[i].x < minx)
      minx = T.vertex[i].x;
    if(T.vertex[i].x > maxx)
      maxx = T.vertex[i].x;
    }

  if( (minz <= 0 && maxz >= 0) && (miny <= 0 && maxy >= 0) && maxx >= 0)
   return 1;
  else
   return 0;
}

int Check_point(float xint, TRIANGLE T)
/* Checks to see if point is inside the triangle */
{
  float angle = 0;
  float v1[4], v2[4], v3[4];

  if(xint >= 0)
    {
    angle = 0.0;
    v1[1] = T.vertex[0].x - xint;
    v1[2] = T.vertex[0].y - 0;
    v1[3] = T.vertex[0].z - 0;
    normalize(v1);

    v2[1] = T.vertex[1].x - xint;
    v2[2] = T.vertex[1].y - 0;
    v2[3] = T.vertex[1].z - 0;
    normalize(v2);

    v3[1] = T.vertex[2].x - xint;
    v3[2] = T.vertex[2].y - 0;
    v3[3] = T.vertex[2].z - 0;
    normalize(v3);

    angle += acos(dot_product(v1, v2));
    angle += acos(dot_product(v2, v3));
    angle += acos((double)dot_product(v3, v1));

    if(fabs(angle - 2*PI) < 0.0001)
      return 1;
    else
      return 0;
    }
  else
    return 0;
}

void Intersect_tri(TRIANGLE T, int organ_id, XP_ARRAY *int_points)
{
  double A,B,C,D;
  double xint, min_error;
  int count, count2;
  int i, j;
  int flag;
  double v1[3], v2[3], v3[3];

  v1[0] = T.vertex[0].x; v1[1] = T.vertex[0].y; v1[2] = T.vertex[0].z;
  v2[0] = T.vertex[1].x; v2[1] = T.vertex[1].y; v2[2] = T.vertex[1].z;
  v3[0] = T.vertex[2].x; v3[1] = T.vertex[2].y; v3[2] = T.vertex[2].z;

  Plane_eqn(v1, v2, v3, &A, &B, &C, &D);
  if(A != 0.0)
    xint = -D/A;

  if(xint <= 0)
    {
    xint = 0.0;
    for(i = 0; i < 3; i++)
      xint += T.vertex[i].x;
    xint = xint / 3.0;
    }

  if(Check_point(xint, T))
    {
    if(int_points->length == 0)
      {
      int_points->length = 1;
      int_points->xp[0].x = xint;
      int_points->xp[0].organ_id = organ_id;
      }
    else
      {
      count = 0;
      flag = 1;
      while(count < int_points->length && flag)
        if(xint > int_points->xp[count].x)
          count++;
        else
          flag = 0;

      for(count2 = int_points->length; count2 > count; count2--)
        {
        int_points->xp[count2].x = int_points->xp[count2-1].x;
        int_points->xp[count2].organ_id = int_points->xp[count2-1].organ_id;
        }

      int_points->xp[count].x = xint;
      int_points->xp[count].organ_id = organ_id;
      int_points->length++;
      }
    }
}

void Find_Intersections_tri(TRI_MODEL *tri_model, int organ_id, float line_origin[3],
float line_vector[3], XP_ARRAY *int_points)
{
  typedef float FLOAT;

  const FLOAT magnitude = sqrt(
	  line_vector[0]*line_vector[0] +
	  line_vector[1]*line_vector[1] +
	  line_vector[2]*line_vector[2]);
  // const float PHI = 180.0f / PI * acos(line_vector[2] / magnitude);
  const FLOAT PHIRad = acos(line_vector[2] / magnitude) - PI / 2.0f;
  // const float THETA = Calc_Angle(line_vector[0], line_vector[1]);
  const FLOAT THETARad = Calc_Angle(line_vector[0], line_vector[1]) * PI / 180.0f;

  // speed optimisation
  const FLOAT cosPhi = cos(PHIRad), sinPhi = sin(PHIRad);
  const FLOAT cosTheta = cos(THETARad), sinTheta = sin(THETARad);

  for(TRIANGLE *tri = tri_model->tris, *triEnd = tri_model->tris + tri_model->num_tris;
      tri != triEnd; tri++) {
    TRIANGLE T;
    for(int j = 0; j < 3; j++) {
	  // double p[3];
      POINT &p = T.vertex[j];
      // point<FLOAT> p;

	  // Translate
      p = tri->vertex[j] - line_origin;

      // Rotate

	  // Rotate_Z(THETA, p);
	  Rotate_Z(cosTheta, sinTheta, p);
	  // Rotate_Y(-(90.0f - PHI), p);
	  Rotate_Y(cosPhi, sinPhi, p);

      // T.vertex[j] = p;
    }

    if(Test_extents_tri(T))
      Intersect_tri(T, organ_id, int_points);
  }
}

int Check_Y_Boundary2(TRI_MODEL slice_tmodel, int x, int y, int z, int intensity)
{
  int j, i, flag;
  float line_origin[3], line_vector[3];
  XP_ARRAY int_points;
  int count_even = 0, count_odd = 0;

  line_origin[0] = (x-xoffset)*pixel_width*10 / subvoxel_index;
  line_origin[1] = (y-yoffset)*pixel_width*10 / subvoxel_index;
  line_origin[2] = (z-zoffset)*slice_width*10 / subvoxel_index;

  line_vector[0] = 1;
  line_vector[1] = 1;
  line_vector[2] = 0;

  int_points.length = 0;
  Find_Intersections_tri(&slice_tmodel, intensity, line_origin, line_vector, &int_points);

  if(int_points.length % 2 == 0)
    count_even++;
  else
    count_odd++;

  line_vector[0] = -1;
  line_vector[1] = -1;
  line_vector[2] = 0;

  int_points.length = 0;
  Find_Intersections_tri(&slice_tmodel, intensity, line_origin, line_vector, &int_points);

  if(int_points.length % 2 == 0)
    count_even++;
  else
    count_odd++;

  if(count_even == count_odd) //Break the tie
    {
    line_vector[0] = 1;
    line_vector[1] = 0;
    line_vector[2] = 0;

    int_points.length = 0;
    Find_Intersections_tri(&slice_tmodel, intensity, line_origin, line_vector, &int_points);

    if(int_points.length % 2 == 0)
      count_even++;
    else
      count_odd++;
    }
/*
  for(j = 0; j < 360; j += 90)
    {
    line_vector[0] = cos(j*PI/180.0);
    line_vector[1] = sin(j*PI/180.0);
    line_vector[2] = 0;

    int_points.length = 0;
    Find_Intersections_tri(&slice_tmodel, intensity, line_origin, line_vector, &int_points);

    if(int_points.length % 2 == 0)
      count_even++;
    else
      count_odd++;
    }
*/

  if(count_odd > count_even)
    return 1;
  else
    return 0;
}

void Fill_tri(TRI_MODEL slice_tmodel, int intensity, float *out, int z)
{
  unsigned long index;
  int edge[200];
  int edge_counter;
  // int counter = 0;
  int i;

  std::fill_n(edge, 200, -1);  // is this really required?

  for (int y = 0; y < tydim; ++y)
  {
    edge_counter = 0;
    float *yslice = out + (y * txdim);
    for (int x = 0; x < txdim; ++x)
    {
      if (int(yslice[x]) == intensity)
        edge[edge_counter++] = x;
    }
    for(int e = edge[--edge_counter], ePrev; edge_counter > 0; --edge_counter)
    {
      ePrev = edge[edge_counter - 1];
      if (e != -1
       && ePrev != -1
       && e - ePrev > 1
	   && Check_Y_Boundary2(slice_tmodel, (ePrev + e) / 2, y, z, intensity))
      {
        G_line(ePrev, tydim - y, e, tydim - y, intensity, out);
      }
	  e = ePrev;
    }
  }
}
