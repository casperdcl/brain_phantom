#include "nurbs.h"
#include "constants.h"
#include "global_vars.h"

#include <malloc.h>
#include <cstdio>
#include <cstdlib>  // exit()
#include <cmath>  // fabs()
#include <memory>  // memcpy()
#include <cassert>

#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
		fprintf(stderr,"Numerical Recipes run-time error...\n");
		fprintf(stderr,"%s\n",error_text);
		fprintf(stderr,"...now exiting to system...\n");
		exit(1);
}

listelement *list_vector(long nl,long nh)
{
		listelement *v = new listelement[nh-nl+1+NR_END];
		if (!v) nrerror("allocation failure in listelement vector()");
		return v-nl+NR_END;
}

void free_list_vector(listelement *v,long nl,long nh)
{
		free((FREE_ARG) (v+nl-NR_END));
		nh=nh;
}

float ***f3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;

  t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
  if (!t) 
	{
	printf("allocation failure 1 in f3tensor()");
	exit(1);
	}
  t +=NR_END;
  t -=nrl;

  t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
  if (!t[nrl]) 
	{
	printf("allocation failure 2 in f3tensor()");
	exit(1);
	}
  t[nrl] +=NR_END;
  t[nrl] -=ncl;

  t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*
  sizeof(float)));
  if (!t[nrl][ncl]) 
	{
	printf("allocation failure 3 in f3tensor()");
	exit(1);
	}
  t[nrl][ncl] +=NR_END;
  t[nrl][ncl] -=ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
	t[i]=t[i-1]+ncol;
	t[i][ncl]=t[i-1][ncl]+ncol*ndep;
	for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  return t;
}

  
void free_f3tensor(float ***t,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
  nrh=nrh;
  nch=nch;
  ndh=ndh;
}

HPOINT *hp_vector(long nl, long nh)
{
  HPOINT *v = new HPOINT[nh-nl+1+NR_END];
  if(!v)
	{
	printf("\nallocation error in hp_vector");
	exit(1);
	}
  return v - nl + NR_END;
}
	 
void free_hpvector(HPOINT *v, long nl, long nh)  
{
  free((FREE_ARG) (v + nl-NR_END));
}
  
float *vector(long nl,long nh)
/* allocate & initialize a float vector with subscript range v[nl..nh] */
{       
		float *v;
		long i;
		v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
		if (!v) nrerror("allocation failure in vector()");
	for (i=nl; i<=nh; i++)
		  v[i]=0.0;
		return v-nl+NR_END;
}
		
void free_vector(float *v,long nl,long nh)
/* free a float vector allocated with vector() */
{
		free((FREE_ARG) (v+nl-NR_END));
		nh=nh;
}

float **matrix(long nrl,long nrh,long ncl,long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
		long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
		float **m;
		 
		/* allocate pointers to rows */
		m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
		if (!m) nrerror("allocation failure 1 in matrix()");
		m += NR_END;
		m -= nrl;
		 
		/* allocate rows and set pointers to them */
		m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
		if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
		m[nrl] += NR_END;
		m[nrl] -= ncl;
		for (i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
		/* return pointer to array of pointers to rows */
		return m;
}
 
void free_matrix(float **m,long nrl,long nrh,long ncl,long nch)
/* free a float matrix allocated by matrix() */
{
		free((FREE_ARG) (m[nrl]+ncl-NR_END));
		free((FREE_ARG) (m+nrl-NR_END));
		nch=nch;
		nrh=nrh;
}

/* Allocation for structures */
TRIANGLE *tri_vector(long nl,long nh)
{
		TRIANGLE *v = new TRIANGLE[nh-nl+1+NR_END];
		if (!v)   
		  {
		  printf("\nallocation failure in tri_vector()");
		  exit(1);
		  } 
	 
		return v-nl+NR_END;
}

void free_tri_vector(TRIANGLE *v,long nl,long nh)
{
		free((FREE_ARG) (v+nl-NR_END));
		nh=nh;
}

/*----------NURBS.C-------------*/
#define TINY 1.0e-20;
#define NR_END 1
#define FREE_ARG char*

BEZIER_PATCH *bp_vector(long nl, long nh)
{
  BEZIER_PATCH *v = new BEZIER_PATCH[nh-nl+1+NR_END];
  if(!v)   
	{
	printf("\nallocation error in bp_vector");
	exit(1);
	}
  return v - nl + NR_END;
}

void free_bpvector(BEZIER_PATCH *v, long nl, long nh)
{
  free((FREE_ARG) (v + nl-NR_END));
}

int *ivector(long nl, long nh)
{
  int *v;

  v = (int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
  if(!v)   
	{
	printf("\nallocation error in ivector");
	exit(1);
	}
  return v - nl + NR_END;
}

void free_ivector(int *v, long nl, long nh)
{
  free((FREE_ARG) (v + nl-NR_END));
}



CPOINT *cp_vector(long nl, long nh)
{
  CPOINT *v = new CPOINT[nh-nl+1+NR_END];
  if(!v) 
	{
	printf("\nallocation error in cp_vector");
	exit(1);
	}
  return v - nl + NR_END;
}

VPOINT *vp_vector(long nl, long nh)
{
  VPOINT *v = new VPOINT[nh-nl+1+NR_END];
  if(!v)
	{
	printf("\nallocation error in vp_vector");
	exit(1);
	}
  return v - nl + NR_END;
}

POINT *p_vector(long nl, long nh)
{
  POINT *v = new POINT[nh-nl+1+NR_END];
  if(!v) 
	{
	printf("\nallocation error in p_vector");
	exit(1);
	}
  return v - nl + NR_END;
}

 
void free_cpvector(CPOINT *v, long nl, long nh)
{
  free((FREE_ARG) (v + nl-NR_END));
}  

void free_pvector(POINT *v, long nl, long nh)
{
  free((FREE_ARG) (v + nl-NR_END));
}  


POINT **p_matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  typedef POINT * Point_;
  POINT **m = new Point_[nrow+NR_END];

  if(!m)
	{
	printf("/n allocation error in p_matrix");
	exit(1);
	}
  m+=NR_END;
  m-=nrl;

  m[nrl]= new POINT[nrow*ncol+NR_END];
  if(!m[nrl])
	{
	printf("/n allocation error in p_matrix");
	exit(1);
	}
  m[nrl] += NR_END;
  m[nrl] -= ncl;
 
  for(i = nrl + 1; i<=nrh; i++) m[i] = m[i-1] + ncol;
  return m;
}


POINT ***p_3d(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  POINT ***t;

  t=(POINT ***) malloc((size_t)((nrow+NR_END)*sizeof(POINT**)));
  if (!t) 
	{
	printf("allocation failure 1 in f3tensor()");
	exit(1);
	}
  t +=NR_END;
  t -=nrl;

  t[nrl]=(POINT **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(POINT*)));
  if (!t[nrl]) 
	{
	printf("allocation failure 2 in f3tensor()");
	exit(1);
	}
  t[nrl] +=NR_END;
  t[nrl] -=ncl;

  t[nrl][ncl]=(POINT *) malloc((size_t)((nrow*ncol*ndep+NR_END)*
  sizeof(POINT)));
  if (!t[nrl][ncl]) 
	{
	printf("allocation failure 3 in f3tensor()");
	exit(1);
	}
  t[nrl][ncl] +=NR_END;
  t[nrl][ncl] -=ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
	t[i]=t[i-1]+ncol;
	t[i][ncl]=t[i-1][ncl]+ncol*ndep;
	for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  return t;
}

  
void free_p_3d(POINT ***t,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
  nrh=nrh;
  nch=nch;
  ndh=ndh;
}

CPOINT ***c_3d(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  CPOINT ***t;

  t=(CPOINT ***) malloc((size_t)((nrow+NR_END)*sizeof(CPOINT**)));
  if (!t) 
	{
	printf("allocation failure 1 in f3tensor()");
	exit(1);
	}
  t +=NR_END;
  t -=nrl;

  t[nrl]=(CPOINT **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(CPOINT*)));
  if (!t[nrl]) 
	{
	printf("allocation failure 2 in f3tensor()");
	exit(1);
	}
  t[nrl] +=NR_END;
  t[nrl] -=ncl;

  t[nrl][ncl]=(CPOINT *) malloc((size_t)((nrow*ncol*ndep+NR_END)*
  sizeof(CPOINT)));
  if (!t[nrl][ncl]) 
	{
	printf("allocation failure 3 in f3tensor()");
	exit(1);
	}
  t[nrl][ncl] +=NR_END;
  t[nrl][ncl] -=ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
	t[i]=t[i-1]+ncol;
	t[i][ncl]=t[i-1][ncl]+ncol*ndep;
	for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  return t;
}

  
void free_c_3d(CPOINT ***t,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
  nrh=nrh;
  nch=nch;
  ndh=ndh;
}


CURVE *curve_vector(long nl, long nh)
{
  CURVE *v;

  v = (CURVE *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(CURVE)));
  if(!v) 
	{
	printf("\nallocation error in curve_vector");
	exit(1);
	}
  return v - nl + NR_END;
}


void free_curve_vector(CURVE *v, long nl, long nh)
{
  free((FREE_ARG) (v + nl-NR_END));
  nh = nh;
}


CURVE **curve_matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  CURVE **m;

  m=(CURVE **) malloc((size_t)((nrow+NR_END)*sizeof(CURVE*)));

  if(!m)
	{
	printf("/n allocation error in curve_matrix");
	exit(1);
	}
  m+=NR_END;
  m-=nrl;

 
  m[nrl]=(CURVE *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(CURVE)));
  if(!m[nrl])
	{
	printf("/n allocation error in curve_matrix");
	exit(1);
	}
  m[nrl] += NR_END;
  m[nrl] -= ncl;
 
  for(i = nrl + 1; i<=nrh; i++) m[i] = m[i-1] + ncol;
  return m; 
}

void free_curve_matrix(CURVE **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}
	

CURVE ***curve_3d(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  CURVE ***t;

  t=(CURVE ***) malloc((size_t)((nrow+NR_END)*sizeof(CURVE**)));
  if (!t) 
	{
	printf("allocation failure 1 in f3tensor()");
	exit(1);
	}
  t +=NR_END;
  t -=nrl;

  t[nrl]=(CURVE **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(CURVE*)));
  if (!t[nrl]) 
	{
	printf("allocation failure 2 in f3tensor()");
	exit(1);
	}
  t[nrl] +=NR_END;
  t[nrl] -=ncl;

  t[nrl][ncl]=(CURVE *) malloc((size_t)((nrow*ncol*ndep+NR_END)*
  sizeof(CURVE)));
  if (!t[nrl][ncl]) 
	{
	printf("allocation failure 3 in f3tensor()");
	exit(1);
	}
  t[nrl][ncl] +=NR_END;
  t[nrl][ncl] -=ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
	t[i]=t[i-1]+ncol;
	t[i][ncl]=t[i-1][ncl]+ncol*ndep;
	for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  return t;
}

  
void free_curve_3d(CURVE ***t,long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
  free((FREE_ARG) (t[nrl]+ncl-NR_END));
  free((FREE_ARG) (t+nrl-NR_END));
  nrh=nrh;
  nch=nch;
  ndh=ndh;
}
  

void free_pmatrix(POINT **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

CPOINT **cp_matrix(long nrl, long nrh, long ncl, long nch)
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  CPOINT **m;
 
  m=(CPOINT **) malloc((size_t)((nrow+NR_END)*sizeof(CPOINT*)));

  if(!m)
	{
	printf("/n allocation error in cp_matrix");
	exit(1);
	}
  m+=NR_END;
  m-=nrl;

  m[nrl]=(CPOINT *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(CPOINT)));
  if(!m[nrl])
	{
	printf("/n allocation error in cp_matrix");
	exit(1);
	}
  m[nrl] += NR_END;
  m[nrl] -= ncl;
 
  for(i = nrl + 1; i<=nrh; i++) m[i] = m[i-1] + ncol;
  return m;
}

void free_cpmatrix(CPOINT **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-NR_END));
  free((FREE_ARG) (m+nrl-NR_END));
}

void ludcmp(float **a, int n, int *indx)
{
	int i,imax,j,k;
	float big,dum,sum,temp;
	float *vv;

	vv=vector(1,n);
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) 

				  {
				   printf("Singular matrix in routine LUDCMP");
				   exit(1);
				  }
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
}

#undef TINY

void lubksb(float **a, int n,int *indx, float b[])
{
	int i,ii=0,ip,j;
	float sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void MakeCurve(CURVE *C, INDEX n, INDEX m, DEGREE p)
{
  int i;

  C->pol.n = n;
  C->pol.Pw = cp_vector(0, n);
  for(i = 0; i <= n; i++)
	{
	C->pol.Pw[i].x = 0.0;
	C->pol.Pw[i].y = 0.0;
	C->pol.Pw[i].z = 0.0;
	C->pol.Pw[i].w = 0.0;
	}
  C->p = p;
  C->knt.m = m;
  C->knt.U = vector(0, m);
  for(i = 0; i <= m; i++)
	C->knt.U[i] = 0.0;

  C->uk = vector(0, n-1);
  for(i = 0; i < n; i++)
	C->uk[i] = 0.0;
}

void BasisFuns(INDEX i, float u, DEGREE p, KNOTVECTOR knot_v, float *N)
{
  INDEX j, r;
  float *left, *right;
  float saved, temp;

  left = vector(1, p);
  right = vector(1, p);
  N[0] = 1.0;
  

  for(j = 1; j<=p; j++)
	{  
	left[j] = u-knot_v.U[i+1-j];
	right[j] = knot_v.U[i+j] - u;
	saved = 0.0;
	for(r = 0; r<j; r++)
	  {
	  temp = N[r]/(right[r+1]+left[j-r]);
	  N[r] = saved+right[r+1]*temp; 
	  saved = left[j-r]*temp;
	  }
	N[j] = saved;
   }
 free_vector(left, 1, p);
 free_vector(right, 1, p);
}

int FindSpan(INDEX n,DEGREE p,float u, KNOTVECTOR knot_v)
{
  int low, high, mid;

  if( (int)(u) == 1) 
	return(knot_v.m - p - 1); 

  low = 0; high = n+1; 
  mid = (low+high) / 2;

  while (u < knot_v.U[mid] || u >= knot_v.U[mid+1])
	{
	if(u < knot_v.U[mid]) 
	  high = mid;
	else         
	  low = mid;
	mid = (low+high) / 2;
	}
  return(mid);
}


POINT CurvePoint(INDEX n, DEGREE p, KNOTVECTOR knot_v, CPOINT *Pw, float u)
{
  CPOINT Cw;
  POINT C;
  float *N;
  int i,j, span;


  N = vector(0, p);
  for(i = 0; i<=p; i++)
	N[i] = 0;
	
  span = FindSpan(n,p,u,knot_v);
  BasisFuns(span,u,p,knot_v,N);

  Cw.x = 0.0; Cw.y = 0.0; Cw.z = 0.0; Cw.w = 0.0;

  for(j = 0; j<=p; j++)
	{ 
	Cw.x = Cw.x + N[j] * Pw[span-p+j].x;
	Cw.y = Cw.y + N[j] * Pw[span-p+j].y;
	Cw.z = Cw.z + N[j] * Pw[span-p+j].z;
	Cw.w = Cw.w + N[j] * Pw[span-p+j].w;
	}

  if(Cw.w != 0)
	{
	C.x = Cw.x / Cw.w;
	C.y = Cw.y / Cw.w;
	C.z = Cw.z / Cw.w;
	}
  else
	{
	C.x = Cw.x;
	C.y = Cw.y;
	C.z = Cw.z;
	}

  free_vector(N, 0, p);
  return C;
}

POINT SurfacePoint(INDEX n, DEGREE p, KNOTVECTOR u_knot, INDEX m, DEGREE q, KNOTVECTOR v_knot, CPOINT **Pw, float u, float v)
{
  CPOINT Sw, *temp;
  POINT S;
  float *Nu, *Nv;
  int l,k, uspan, vspan;

  temp = cp_vector(0, q);

  Nu = vector(0, p);
  Nv = vector(0, q);

  for(l = 0; l<=p; l++)
	Nu[l] = 0;

  for(l = 0; l<=q; l++)
	Nv[l] = 0;

  uspan = FindSpan(n,p,u,u_knot);
  BasisFuns(uspan, u,p,u_knot,Nu);

  vspan = FindSpan(m,q,v,v_knot);
  BasisFuns(vspan, v,q,v_knot,Nv);

  for(l = 0; l<=q; l++)
	{ 
	temp[l].x = 0.0;
	temp[l].y = 0.0;
	temp[l].z = 0.0;
	temp[l].w = 0.0;

	for(k = 0; k<=p; k++) 
	  {
	  temp[l].x = temp[l].x + Nu[k] * Pw[uspan-p+k][vspan-q+l].x;
	  temp[l].y = temp[l].y + Nu[k] * Pw[uspan-p+k][vspan-q+l].y;
	  temp[l].z = temp[l].z + Nu[k] * Pw[uspan-p+k][vspan-q+l].z;  
	  temp[l].w = temp[l].w + Nu[k] * Pw[uspan-p+k][vspan-q+l].w; 
	  } 
	}


  Sw.x = 0.0; Sw.y = 0.0; Sw.z = 0.0; Sw.w = 0.0;

  for(l = 0; l<=q; l++)
	{
	Sw.x = Sw.x + Nv[l] * temp[l].x;
	Sw.y = Sw.y + Nv[l] * temp[l].y;
	Sw.z = Sw.z + Nv[l] * temp[l].z;
	Sw.w = Sw.w + Nv[l] * temp[l].w;
	}


  if(Sw.w !=0)
	{
	S.x = Sw.x / Sw.w;
	S.y = Sw.y / Sw.w;
	S.z = Sw.z / Sw.w;
	}
  else
	{
	S.x = Sw.x;
	S.y = Sw.y;
	S.z = Sw.z;
	}


  free_cpvector(temp,0, q);
  free_vector(Nu, 0, p);
  free_vector(Nv, 0, q);

  return S;
}


float distance3d(POINT p1, POINT p2)
{
  float tempx, tempy, tempz;

  tempx = p1.x - p2.x;
  tempy = p1.y - p2.y;
  tempz = p1.z - p2.z;
	
  return( sqrt(tempx * tempx + tempy * tempy + tempz * tempz));
}


void Calc_UniformKnotVector(INDEX n, DEGREE p, KNOTVECTOR knot_v, float *uk)
{
  INDEX i,j, m;
  float d = 0.0;
  float sum = 0.0;
  
  m = n+p+1;
  
  uk[0] = 0; uk[n] = 1;

  for(i = 1; i < n; i++)
	uk[i] = (float)i/(float)n;
  
/* Calculate Knot Vector */
  for(i = 0; i<=p; i++)
	  knot_v.U[i] = 0;
  for(i = m-p; i<= m; i++)
	  knot_v.U[i] = 1;
  for(j = 1; j <= n - p; j++)
	{
	sum = 0.0;
	for(i = j; i <= j + p - 1; i++)
	  sum += uk[i];
	knot_v.U[j + p] = sum / p;
	}
	
  knot_v.m = m;
}

void Calc_KnotVector(INDEX n, DEGREE p, POINT *Qw, KNOTVECTOR knot_v, float *uk)
{
  INDEX i,j, m;
  float d = 0.0;
  float sum = 0.0;
  
  m = n+p+1;
  
  uk[0] = 0; uk[n] = 1;
  for(i = 1; i<= n; i++)
	d += distance3d(Qw[i], Qw[i - 1]);

  if(d == 0)
	{
	for(i = 1; i < n; i++)
	  uk[i] = (float)i/(float)n;
	}

  else
	{
	for(i = 1; i< n; i++)
	  uk[i] = uk[i - 1] + (distance3d(Qw[i], Qw[i - 1]) / d);
	}
  
/* Calculate Knot Vector */
  for(i = 0; i<=p; i++)
	  knot_v.U[i] = 0;
  for(i = m-p; i<= m; i++)
	  knot_v.U[i] = 1;
  for(j = 1; j <= n - p; j++)
	{
	sum = 0.0;
	for(i = j; i <= j + p - 1; i++)
	  sum += uk[i];
	knot_v.U[j + p] = sum / p;
	}
	
  knot_v.m = m;
}


void GlobalCurveInterp(INDEX n, POINT *Qw, int r, DEGREE p, KNOTVECTOR knot_v, float *uk, CPOINT *Pw)
{
  float **mA, **mA2;
  float *rhs, *ipointer;
  int *indx;
  float d = 0.0;
  float sum = 0.0;
  INDEX i, j;
  int span;

  mA = matrix(0, n, 0, n);
  mA2 = matrix(1, n+1, 1, n+1);
  indx = ivector(1, n+1);
  rhs = vector(1, n+1);

/*Initialize mA to zero*/
  for(i = 0; i<=n; i++)
	for(j = 0; j<=n; j++)
	  mA[i][j] = 0.0;

  for(i = 0; i<=n; i++)
   { 
   span = FindSpan(n, p, uk[i], knot_v);
   ipointer = &mA[i][span - p];
   BasisFuns(span, uk[i], p, knot_v, ipointer);
   }  

  for(i = 0; i<=n; i++)
   {
   for(j = 0; j<=n; j++)
	 {
	 mA2[i+1][j+1] = mA[i][j];
	 }
   }

  ludcmp(mA2, n+1, indx);

  for(i = 0; i<r; i++)
	{
	for(j = 0; j<=n; j++)  
	  {
	  if(i == 0)
		rhs[j+1] = Qw[j].x;
	  else if(i == 1)
		rhs[j+1] = Qw[j].y;
	  else if(i == 2)
		rhs[j+1] = Qw[j].z;
	  }
	lubksb(mA2, n+1, indx, rhs); 
	for(j = 0; j<=n; j++)  
	  {
	  if(i == 0)
		Pw[j].x = rhs[j+1];
	  if(i == 1)
		Pw[j].y = rhs[j+1];
	  if(i == 2)
		Pw[j].z = rhs[j+1];
	  } 
	}

  free_ivector(indx, 1, n+1);
  free_vector(rhs, 1, n+1);
  free_matrix(mA, 0, n, 0, n); 
  free_matrix(mA2, 1, n+1, 1, n+1); 
}


void SurfMeshParams(INDEX n, INDEX m, QNET Q, float *uk, float *vl)
{
  INDEX num, k, l, max;
  float total, *cds, d;

  if(m >= n) 
	max = m;
  else
	max = n;

  cds = vector(1, max + 1);
  num = m+1;

  uk[0] = 0.0; uk[n] = 1.0;

  for(k = 1; k<n; k++) uk[k] = 0.0;

  for(l = 0; l<=m; l++)
	{
	total = 0.0;

	for(k = 1; k<=n; k++)
	  {
	  cds[k] = distance3d(Q.Qw[k][l], Q.Qw[k-1][l]);
	  total = total + cds[k];
	  }

	if(total == 0.0) num = num - 1;
	else
	  {
	  d = 0.0; 
	  for(k = 1; k<n; k++)
		{
		d = d + cds[k];
		uk[k] = uk[k] + d / total;
		}
	  }
	}

  if(num == 0) 
	{
	printf("\n Error in SurfMeshParams!!!");
	exit(1);
	}

  for(k = 0; k<n; k++)  
	uk[k] = uk[k] / num; 

  num = n + 1;

  vl[0] = 0.0; vl[m] = 1.0;
  for(l = 1; l<m; l++) vl[l] = 0.0;
  for(k = 0; k<=n; k++)
	{
	total = 0.0;
	for(l = 1; l<=m; l++)
	  {
	  cds[l] = distance3d(Q.Qw[k][l], Q.Qw[k][l-1]);
	  total = total + cds[l];
	  }
	if(total == 0.0) num = num - 1;
	else
	  {
	  d = 0.0; 
	  for(l = 1; l<m; l++)
		{
		d = d + cds[l];
		vl[l] = vl[l] + d / total;
		}
	  }
	}

  if(num == 0) 
	{
	printf("\n Error in SurfMeshParams!!!");
	exit(1);
	}

  for(l = 0; l<m; l++)  
	vl[l] = vl[l] / num;

  free_vector(cds, 1, max + 1);
}
  

void GlobalSurfInterp(INDEX n, INDEX m, QNET Q, DEGREE p, DEGREE q, KNOTVECTOR u_knot, KNOTVECTOR v_knot, CNET P)
{
  INDEX i, j, l, k, max, uknots, vknots;
  float sum;
  QPOINTS QP;
  CNET R;
  float *uk, *vl;
  CPOINT **temp, *Rpw;
  
  temp = cp_matrix(0, n, 0, m);
  uk = vector(0, n+1);
  vl = vector(0, m+1);

  if(m >= n)
	max = m;
  else
	max = n;

  Rpw = cp_vector(0, max);
  QP.Qw = p_vector(0, max);

  for( i = 0; i<=n; i++)
	for(j = 0; j<=m; j++)
	  {
	  temp[i][j].x = 0;
	  temp[i][j].y = 0;
	  temp[i][j].z = 0;
	  temp[i][j].w = 0;
	  }

  R.Pw = temp;
  SurfMeshParams(n, m, Q, uk, vl); 
  uknots = n + p + 1;
  vknots = m + q + 1;

/* Calculate Knot Vector U */
  for(i = 0; i<=p; i++)
	u_knot.U[i] = 0;
  for(i = uknots-p; i<= uknots; i++)
	u_knot.U[i] = 1;
  for(j = 1; j <= n - p; j++)
	{
	sum = 0.0;
	for(i = j; i <= j + p - 1; i++)
	  sum += uk[i];
	u_knot.U[j + p] = sum / p;
	}

/* Calculate Knot Vector V */
  for(i = 0; i<=q; i++)
	v_knot.U[i] = 0;
  for(i = vknots-q; i<= vknots; i++)
	v_knot.U[i] = 1;
  for(j = 1; j <= m - q; j++)
	{
	sum = 0.0;
	for(i = j; i <= j + q - 1; i++)
	  sum += vl[i];
	v_knot.U[j + q] = sum / q;
	}

  for(l = 0; l<=m; l++)
   {
	 for(k = 0; k<=n; k++)
	   {
	   QP.Qw[k].x = Q.Qw[k][l].x;        
	   QP.Qw[k].y = Q.Qw[k][l].y;        
	   QP.Qw[k].z = Q.Qw[k][l].z;
	   }
	 GlobalCurveInterp(n, QP.Qw, 3, p, u_knot, uk, Rpw); 
	 for(i = 0; i<=n; i++)
	   {
	   R.Pw[i][l].x = Rpw[i].x;
	   R.Pw[i][l].y = Rpw[i].y;
	   R.Pw[i][l].z = Rpw[i].z;
	   R.Pw[i][l].w = Rpw[i].w;
	   }
   }

  for(i = 0; i<=n; i++)
   {
	 for(l = 0; l<=m; l++)
	   {
	   QP.Qw[l].x = R.Pw[i][l].x;        
	   QP.Qw[l].y = R.Pw[i][l].y;        
	   QP.Qw[l].z = R.Pw[i][l].z;
	   }
	 GlobalCurveInterp(m, QP.Qw, 3, q, v_knot, vl, Rpw);
	 for(l = 0; l<=m; l++)
	   {
	   P.Pw[i][l].x = Rpw[l].x;
	   P.Pw[i][l].y = Rpw[l].y;
	   P.Pw[i][l].z = Rpw[l].z;
	   P.Pw[i][l].w = Rpw[l].w; 
	   }
   }

  free_cpmatrix(temp,0, n, 0, m);

  free_cpvector(Rpw, 0, max);
  free_pvector(QP.Qw, 0, max);
}


/*-----------------------------------Bezier Patch routines-------------------------------------------*/
/* The following routines are used in converting a NURBS surface into Bezier patches                 */
/*---------------------------------------------------------------------------------------------------*/
int get_breakpoint(int len_kU, float *kU, float w)
{
  register int i;
  i = 0;
  while ((i < len_kU) && (kU[i] <= w))
	i++;
  return (i-1);
}

static struct knotmultCnt_s
{
   struct knotmult_s
	 {
	 int mu;                    /* multiplicity */
	 float val;                 /* value */
	 } *umult, *vmult;
  int numU, numV;               /* # distinct knots */
} knots;

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
void refine_patch(patch *src, patch *dest)
/*  This routine refines the Bezier patches created for a NURBS surface */
{
  register int i, j, s, i2;
  int k, last, r, li, delta;
  float *kU, *kW, *kV;
  CPOINT C[4], *D;
  float lenU, lenV, omega;

  /* Refine along U-direction first */
  if (dest->numU > src->numU)   /* there are new breakpoints */
	{
	k = src->ordU;
	kU = src->kU;
	kW = dest->kU;
	lenU = src->numU + src->ordU; /* #knots in U-direction */
	  
	for (i = 0; i < src->numV; i++)  /* for each V-vertex row do */
	  for (j = 0; j < dest->numU; j++)
		{                       /* for each w[j], calculate the new */
								/* vertex  */
		delta = get_breakpoint(lenU, kU, kW[j]);

		for (s = 0; s <= MIN(k - 1, delta); s++)
		  C[s] = src->points[i][delta - s];     /* Initialize C-array */

		/* Calculate the new vertex W[j] for this value of kW[j] */
		for (r = k; r > 1; r--)
		  {
		  li = delta; 
		  last = MIN(r - 2, delta);
		  for (s = 0; s < last; s++)
			{
			omega = (kW[j + r - 1] - kU[li])/ (kU[li + r - 1] - kU[li]);
			C[s].x = omega*C[s].x + (1 - omega)*C[s+1].x;
			C[s].y = omega*C[s].y + (1 - omega)*C[s+1].y;
			C[s].z = omega*C[s].z + (1 - omega)*C[s+1].z;
			C[s].w = omega*C[s].w + (1 - omega)*C[s+1].w;
			li--; 
			}
		  omega = (kW[j + r - 1] - kU[li])/ (kU[li + r - 1] - kU[li]);
		  if (last < (r - 2))
			{
			C[s].x = omega*C[s].x;
			C[s].y = omega*C[s].y;
			C[s].z = omega*C[s].z;
			C[s].w = omega*C[s].w;
			}
		  else
			{
			C[s].x = omega*C[s].x + (1 - omega)*C[s+1].x;
			C[s].y = omega*C[s].y + (1 - omega)*C[s+1].y;
			C[s].z = omega*C[s].z + (1 - omega)*C[s+1].z;
			C[s].w = omega*C[s].w + (1 - omega)*C[s+1].w;
			}
		  }  
		/* Now, the value of the new vertex is available in C[0] */
		dest->points[i][j] = C[0];
		}
	}
  else                          /* no new breakpoints in U-direction */
	for (i = 0; i < src->numV; i++)
	  memcpy(dest->points[i], src->points[i], src->numU*sizeof(CPOINT));
	  
  /* Now perform refinement along the V-direction */
  if (dest->numV > src->numV)   /* there are new breakpoints */
	{
	k = src->ordV;
	kV = src->kV;
	kW = dest->kV;
	lenV = src->numV + src->ordV; /* #knots in V-direction */
			 
	D = (CPOINT *) malloc(src->numV*(sizeof(CPOINT)));
			 
	for (i = 0; i < dest->numU; i++)  /* for each U-vertex col do */
	  {
	  for (i2 = 0; i2 < src->numV; i2++) /* Make a copy of the */
		D[i2] = dest->points[i2][i];     /* i'th column */
			 
	  for (j = 0; j < dest->numV; j++)
		{               /* for each w[j], calculate the new */
						/* vertex  */
		delta = get_breakpoint(lenV, kV, kW[j]);
	 
		for (s = 0; s <= MIN(k - 1, delta); s++)
	 C[s] = D[delta - s];/* Initialize C-array */
	  
	  /* Calculate the new vertex W[j] for this value of kW[j] */
	  for (r = k; r > 1; r--)
		{
		li = delta;
		last = MIN(r - 2, delta);
		for (s = 0; s < last; s++)
		  {
		  omega = (kW[j + r - 1] - kV[li])/ (kV[li + r - 1] - kV[li]);
		  C[s].x = omega*C[s].x + (1 - omega)*C[s+1].x;
		  C[s].y = omega*C[s].y + (1 - omega)*C[s+1].y;
			   C[s].z = omega*C[s].z + (1 - omega)*C[s+1].z;
		  C[s].w = omega*C[s].w + (1 - omega)*C[s+1].w;
		  li--;
		  }
		omega = (kW[j + r - 1] - kV[li])/ (kV[li + r - 1] - kV[li]);
		if (last < (r - 2))
		  {
		  C[s].x = omega*C[s].x;
		  C[s].y = omega*C[s].y;
			   C[s].z = omega*C[s].z;
			   C[s].w = omega*C[s].w;
		  }
		else
		  {
		  C[s].x = omega*C[s].x + (1 - omega)*C[s+1].x;
		  C[s].y = omega*C[s].y + (1 - omega)*C[s+1].y;
		  C[s].z = omega*C[s].z + (1 - omega)*C[s+1].z;
		  C[s].w = omega*C[s].w + (1 - omega)*C[s+1].w;
		  }
		}
	  /* Now, the value of the new vertex is available in C[0] */
	  dest->points[j][i] = C[0];
	  }
	}
	free(D);
  }
  free(knots.umult);
  free(knots.vmult);
}

#define KNOTEPS		1e-15
#define EPSEQ(a,b)	(fabs ((float)(a-b)) < KNOTEPS)
int get_knot_multiplicities(patch *src)
/* This function checks the multiplicity for each knot of a NURBS surface */
/* To convert a NURBS surface into Bezier patches, each knot must have a multiplicity of 4*/
{
  register int i;
  int curr, lenU, lenV, ncurr;
  float u, v;
		  
  i = curr = 0;
  lenU =  src->ordU + src->numU;
  knots.umult = new knotmultCnt_s::knotmult_s[lenU];

  while (i < lenU)      /* traverse all the U-knots */
	{
	u = src->kU[i];
	for (ncurr = 0; (i < lenU) && EPSEQ(src->kU[i], u); i++, ncurr++)
	  ;
	knots.umult[curr].val = u;
	if(ncurr > src->ordU) ncurr = src->ordU;
	  knots.umult[curr].mu = ncurr;
	curr++;
	}   
  knots.numU = curr-1;

  i = 0;
  curr = 0;
  lenV =  src->ordV + src->numV;
  knots.vmult = new knotmultCnt_s::knotmult_s[lenV];

  while (i < lenV)      /* traverse all the V-knots */
	{
	v = src->kV[i];
	for (ncurr = 0; (i < lenV) && EPSEQ(src->kV[i], v); i++, ncurr++)
	  ;
	knots.vmult[curr].val = v;
	if(ncurr > src->ordV) ncurr = src->ordV;
	  knots.vmult[curr].mu = ncurr;
	curr++;
	}
  knots.numV = curr-1;
  return 0;
}

int alloc_patch(patch *p)
/* Allocate space for storing knot arrays and the vertices */
/* constituting the patch */
{
  int i;
 
  if (((p->kU = (float *) malloc((p->numU + p->ordU)*sizeof(float)))
	   ==  NULL) ||
	  ((p->kV = (float *) malloc((p->numV + p->ordV)*sizeof(float)))
	   ==  NULL) ||
	  ((p->points = (CPOINT **) malloc(p->numV*sizeof(CPOINT *))) ==
	   NULL))
	{
	perror("malloc");
	exit(1);
	}

  for (i = 0; i < p->numV; i++)
	 if ((p->points[i] = (CPOINT *) malloc(p->numU*sizeof(CPOINT))) == NULL)
	   {
	   perror("malloc");
	   exit(1);
	   }
   return 0;
}

int free_patch(patch *p)
/* Free the space allocated for this patch since it is no longer */
/* required */
{
  int i;
	   
  free(p->kU);
  free(p->kV);
  for (i = 0; i < p->numV; i++)
   free(p->points[i]);
	 
  free(p->points);   
  return 0; 
}


#define MAXU 1000
#define MAXV 1000
#define MAXT 1000
int insert_multiple_knots(patch *src, patch *dest)
/* Function inserts knots into a NURBS surface: Knots are inserted into a NURBS surface in order to */
/* convert it into Bezier patches */
{
  register int i, j;
  int ucurr, vcurr, kcurr;
  float tmpkU[MAXU], tmpkV[MAXV];

  assert(src->ordU <= MAX_ORDER);
  assert(src->ordV <= MAX_ORDER);
	   
  *dest = *src;
  get_knot_multiplicities(src);
  
  dest->kU = tmpkU;   
  dest->kV = tmpkV;
  
  kcurr = ucurr = i = 0;
 
  while (i + knots.umult[kcurr].mu < dest->ordU)
	{
	for (j = 0; (j < knots.umult[kcurr].mu); j++, i++)
	  dest->kU[ucurr++] = knots.umult[kcurr].val;
	kcurr++;
	}
  while (i < src->numU + 1)
	{
	for (j = 0; j < dest->ordU; j++)
	  dest->kU[ucurr++] = knots.umult[kcurr].val;
  
	i += knots.umult[kcurr].mu;  
	kcurr++;
	}
  while (i < src->numU + src->ordU)
	 dest->kU[ucurr++] = src->kU[i++];

  /* Repeat the same process for V */
  kcurr = vcurr = i = 0;
  
  while (i + knots.vmult[kcurr].mu < dest->ordV)
	{
	for (j = 0; (j < knots.vmult[kcurr].mu); j++, i++)
	  dest->kV[vcurr++] = knots.vmult[kcurr].val;
	kcurr++;
	}
  while (i < src->numV + 1)
	{
	for (j = 0; j < dest->ordV; j++)
	  dest->kV[vcurr++] = knots.vmult[kcurr].val;
  
	i += knots.vmult[kcurr].mu;
	kcurr++;
	}
  while (i < src->numV + src->ordV)
	 dest->kV[vcurr++] = src->kV[i++];

  /* Now allocate space for new patch to hold the knot arrays and the */
  /* vertices */
  dest->numU = ucurr - dest->ordU;
  dest->numV = vcurr - dest->ordV;
  alloc_patch(dest);
  memcpy(dest->kU, tmpkU, ucurr*sizeof(float));
  memcpy(dest->kV, tmpkV, vcurr*sizeof(float));
  return 0;
}

int setup_initial_patch(patch *p, SURFACE *nrb_model)
/* This routine sets up the initial Bezier patch (p) for the NURBS surface (nrb_model) */
{
  register int i, j;
  int lenU, lenV;
  float *tmpkU;
	 
  p->ordU = 3; p->ordV = 3;
  p->ordU++; p->ordV++;
	
  lenU = nrb_model->net.n + 4;
  p->numU  = lenU - p->ordU;
  assert((tmpkU = (float *)malloc(lenU * sizeof(float))) != NULL);
  for (i = 0; i < lenU; i++) /* Read in the U-knots */
	tmpkU[i] = nrb_model->knu.U[i];
  
  lenV = nrb_model->net.m + 4;
  p->numV  = lenV - p->ordV;
  alloc_patch(p);              /* Allocate space for the tables */
  memcpy(p->kU, tmpkU, lenU*sizeof(float));
  free(tmpkU);
  for (i = 0; i < lenV; i++) /* Read in the V-knots */
	 p->kV[i] = nrb_model->knv.U[i];
  
  lenU = nrb_model->net.n;
  lenV = nrb_model->net.m;
 
  assert(lenU == p->numU && lenV == p->numV);
  
  for (i = 0; i < p->numU; i++)   /* read rational vertices */
	 for (j = 0; j < p->numV; j++)
	   {
	   p->points[j][i].x = nrb_model->net.Pw[i][j].x;
	   p->points[j][i].y = nrb_model->net.Pw[i][j].y;
	   p->points[j][i].z = nrb_model->net.Pw[i][j].z;
	   p->points[j][i].w = 1;
	   }
  return 1;
}
  
int create_bezier_patches(patch *p, BEZIER_PATCH *patches)
/* Routine decomposes a NURBS surface (initial patch p) into many Bezier patches (patches) */
/* This is done by inserting knots into the NURBS surface until each knot has a multiplicity of 4 */
{
  register int i, j, k, l, l1, k1;
  int umark, vmark;
  int numpt1, numpt2;
  float uval, vval;
  extern statistics_t pstat;
  float minz, maxz, miny, maxy, minx, maxx;
  float pixel_factor;
  pstat.nbezs = 0;
  
  pixel_factor = subvoxel_index / (slice_width * 10.0);
	 
  uval = p->kU[p->ordU - 1];
  umark = p->ordU - 1;
  while ((umark > 0) && (p->kU[umark] == uval))
	umark--;
  if (p->kU[umark] < uval)   
	umark++;

  vval = p->kV[p->ordV - 1];
  vmark = p->ordV - 1;
  while ((vmark > 0) && (p->kV[vmark] == vval))
	vmark--;
  if (p->kV[vmark] < vval)
	vmark++;
  
  for (l1=0, l = umark; l <= (p->numU - p->ordU); l += p->ordU,l1++)
	{
	for (k1=0, k = vmark; k <= (p->numV - p->ordV); k += p->ordV,k1++)
	  {
	  minz = p->points[k][l].z;
	  maxz = minz;
	  miny = p->points[k][l].y;
	  maxy = miny;
	  minx = p->points[k][l].x;
	  maxx = minx;
	  for (j = 0; j < p->ordU; j++)
		for (i = 0; i < p->ordV; i++)
		  {
		  patches[pstat.nbezs].cntrl_points[j][i][0] = p->points[k+i][l+j].x;
		  patches[pstat.nbezs].cntrl_points[j][i][1] = p->points[k+i][l+j].y;
		  patches[pstat.nbezs].cntrl_points[j][i][2] = p->points[k+i][l+j].z;

		  patches[pstat.nbezs].slice_points[j][i][0] = (p->points[k+i][l+j].x*subvoxel_index) / (10*pixel_width) + xoffset;
		  patches[pstat.nbezs].slice_points[j][i][1] = (p->points[k+i][l+j].y*subvoxel_index) / (10*pixel_width) + yoffset;
		  patches[pstat.nbezs].slice_points[j][i][2] = (p->points[k+i][l+j].z*subvoxel_index) / (10*slice_width) + zoffset;

		  patches[pstat.nbezs].slice_points[j][i][1] = tydim - patches[pstat.nbezs].slice_points[j][i][1];

		  if(p->points[k+i][l+j].z > maxz)
			maxz = p->points[k+i][l+j].z;
		  if(p->points[k+i][l+j].z < minz)
			minz = p->points[k+i][l+j].z;
		  if(p->points[k+i][l+j].y > maxy)
			maxy = p->points[k+i][l+j].y;
		  if(p->points[k+i][l+j].y < miny)
			miny = p->points[k+i][l+j].y;
		  if(p->points[k+i][l+j].x > maxx)
			maxx = p->points[k+i][l+j].x;
		  if(p->points[k+i][l+j].x < minx)
			minx = p->points[k+i][l+j].x;
		  }
  
	   patches[pstat.nbezs].slice_minz = floor(minz*pixel_factor + 0.5);
	   patches[pstat.nbezs].slice_maxz = floor(maxz*pixel_factor + 0.5);

	   patches[pstat.nbezs].maxx = maxx;
	   patches[pstat.nbezs].minx = minx;
	   patches[pstat.nbezs].maxy = maxy;
	   patches[pstat.nbezs].miny = miny;
	   patches[pstat.nbezs].maxz = maxz;
	   patches[pstat.nbezs].minz = minz;    
	   numpt1 = 0; numpt2 = 0;

	   pstat.nbezs++;
	   pstat.tot_trimpts1 += numpt1;
	   pstat.tot_trimpts2 += numpt2;
	  }
	}
  return 0;
}

statistics_t pstat;

void SETUP_BEZIER_MODEL(SURFACE surf, BEZIER_MODEL *bez_model)
{
  bez_model->num_patches = (surf.net.n-3) * (surf.net.m-3);
  bez_model->patches = bp_vector(0, bez_model->num_patches);
}

void SPLINE2BEZ(SURFACE *nrb_model, BEZIER_MODEL *bez_model)
/* This is the main routine for converting a NURBS surface (nrb_model) into a Bezier model (bez_model) */
{
  patch src, dest;

  pstat.count = 0;
  pstat.tot_nbezs = 0;
  pstat.tot_trimpts1 = 0;
  pstat.tot_trimpts2 = 0;

  setup_initial_patch(&src, nrb_model);
  insert_multiple_knots(&src, &dest);
  refine_patch(&src, &dest);
  create_bezier_patches(&dest, bez_model->patches);

  free_patch(&dest);
  free_patch(&src);
}
