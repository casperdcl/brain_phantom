#pragma once

#include "nurbs.h"

void Fill_tri(TRI_MODEL slice_tmodel, int intensity, float *out, int z);
void GET_DYN_PARAMS(char *parfile);
void MoveTriangleModel(TRI_MODEL *tmodel, float x_trans, float y_trans, float z_trans, float x_rot, float y_rot, float z_rot, float *tx, float *ty, float *tz);
void Read_motion_file(char *file, CURVE *C);
void Read_Triangle_ModelVTK(char *filename, TRI_MODEL *tmodel);
void Render_triangle(TRIANGLE tri, int intensity, float *out, int z);
int Test_extents_tri_z(TRIANGLE tri, int z);
