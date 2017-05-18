#pragma once

#include <stdexcept>

typedef long INDEX;
typedef short FLAG;
typedef int INTEGER;
typedef char * STRING;
typedef short DEGREE;
typedef double PARAMETER;

#define NR_END 1
#define SWAP(a,b) itemp=(a); (a) = (b); (b) = itemp;
#define NSTACK 50
#define FREE_ARG char*


template<typename F=float>
struct point
{
  F x;
  F y;
  F z;

  F &operator[](int i) {
  #if !defined(NDEBUG) && (defined(_DEBUG) || defined(DEBUG))
    if(0 > i || i > 2)
      throw std::length_error("Index out of bounds");
  #endif
    return reinterpret_cast<F *>(this)[i];
  }

  const F &operator[](int i) const {
  #if !defined(NDEBUG) && (defined(_DEBUG) || defined(DEBUG))
    if(0 > i || i > 2)
      throw std::length_error("Index out of bounds");
  #endif
    return reinterpret_cast<const F *>(this)[i];
  }

  template<typename FloatArrayType>
  point<F> &operator=(const FloatArrayType &other) {
    x = other[0];
    y = other[1];
    z = other[2];
    return *this;
  }

  template<typename FloatArrayType>
  point<F> operator-(const FloatArrayType &other) const {
    point<F> res;
    res.x = x - other[0];
    res.y = y - other[1];
    res.z = z - other[2];
    return res;
  }

  // cross product
  template<typename FloatArrayType>
  point<F> operator^(const FloatArrayType &other) const {
    point<F> res;
    res.x = y * other[2] - z * other[1];
    res.y = z * other[0] - x * other[2];
    res.z = x * other[1] - y * other[0];
    return res;
  }

  // dot product
  template<typename FloatArrayType>
  F operator*(const FloatArrayType &other) const {
    return x * other[0] + y * other[1] + z * other[2];
  }

  F magSq() const { return x*x + y*y + z*z; }
  F mag() const { return sqrt(magSq()); }

  template<typename FloatArrayType>
  F angle(const FloatArrayType &other) const {
	return acos(*this * other /
      (sqrt(magSq() *
       ( other[0] * other[0]
       + other[1] * other[1]
       + other[2] * other[2]))));
  }
};
typedef point<float> POINT;

typedef struct list_node{
  POINT o_start;
  POINT o_end;

  POINT start;
  POINT end;
  float radius;
  int ID;
  int visible;

  int sZ1, sZ2, sY1, sY2, sX1, sX2;
  int eZ1, eZ2, eY1, eY2, eX1, eX2;
  
  float s_scalez, s_scaley, s_scalex;
  float e_scalez, e_scaley, e_scalex;
} listelement;


typedef struct hpoint
{
  float u;
  float v;
  float p;

  int   id;
  int edge;
  int count;
} HPOINT;

typedef struct bezier_patch
{
  double cntrl_points[4][4][3];
  double slice_points[4][4][3];

  int slice_minz, slice_maxz;
  double minz, maxz, miny, maxy, minx, maxx;
} BEZIER_PATCH;

typedef struct bezier_model
{
  BEZIER_PATCH *patches;
  int num_patches;
} BEZIER_MODEL;


typedef struct qpoints
{
  INDEX n;
  POINT *Qw;
} QPOINTS;

typedef struct cpoint
{
  float  x,
         y,
         z, 
         w;
} CPOINT;

typedef struct vpoint
{
  float  x,
         y,
         z;
  int im1;
  POINT orig;
  int Z1, Z2, Y1, Y2, X1, X2;
  float scalez, scaley, scalex;

  float prev_x, prev_y, prev_z;

  int flag;
} VPOINT;

typedef struct cpolygon
{
  INDEX n;
  CPOINT *Pw;
} CPOLYGON;

typedef struct knotvector
{
    INDEX m;
    float *U;
} KNOTVECTOR;

typedef struct curve
{
  CPOLYGON pol;
  DEGREE p;
  KNOTVECTOR knt;
  float *uk;
} CURVE;

typedef struct cnet
{
  INDEX n,
        m;
  CPOINT **Pw;
} CNET;

typedef struct qnet
{
  INDEX n,
        m;
  POINT **Qw;
} QNET;

typedef struct surface
{
  CNET net;
  DEGREE p,
         q;
  KNOTVECTOR knu, 
             knv;
  int max_z, min_z;
  float xy_span, z_span;
  int flag;
} SURFACE;

typedef struct xpoint
{
  double x;
  int organ_id;
} XPOINT;

typedef struct xp_array
{
  XPOINT xp[100];
  int length;
} XP_ARRAY;

typedef struct triangle
{
  POINT vertex[3];
} TRIANGLE;

typedef struct tri_model
{
  TRIANGLE *tris;
  int num_tris;
  float min_x, max_x, min_y, max_y, min_z, max_z;
} TRI_MODEL;


/*NURBS_BEZ.H*/
#include <stdio.h>
#include <assert.h>

#define MAX_ORDER 20

typedef unsigned int boolean;

/* Declaration of a Nurb surface patch */
typedef struct {
  int             numU, numV;   /* #points in U and V directions */
  int             ordU, ordV;   /* order of the spline in U and V */
  float          *kU, *kV;      /* knot vectors */
                                /* length(kU) = [0...numU - 1 + ordU] */
                                /* length(kV) = [0...numV - 1 + ordV] */
  CPOINT         **points;      /* [0..numV - 1][0..numU -1] array */
} patch;

/* Structure to hold the statistics */
typedef struct stat_s {
  int count;                    /* number of patches */
  int nbezs;                    /* #Bezier patches generated */
  int tot_nbezs;                /* Total #Bez patches for the whole */
  int tot_trimpts1;             /* Max PW Trim Curve Size */
  int tot_trimpts2;             /* Max SPLINE Trim Curve Size */
                                /* file  */
} statistics_t;

void Calc_UniformKnotVector(INDEX n, DEGREE p, KNOTVECTOR knot_v, float *uk);
CPOINT **cp_matrix(long nrl, long nrh, long ncl, long nch);
POINT CurvePoint(INDEX n, DEGREE p, KNOTVECTOR knot_v, CPOINT *Pw, float u);
void free_ivector(int *v, long nl, long nh);
void free_pvector(POINT *v, long nl, long nh);
int *ivector(long nl, long nh);
void MakeCurve(CURVE *C, INDEX n, INDEX m, DEGREE p);
POINT *p_vector(long nl, long nh);
TRIANGLE *tri_vector(long nl,long nh);
TRIANGLE *tri_vector(long nl,long nh);
float *vector(long nl,long nh);
