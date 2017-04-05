/*---------------------------------------------------------------------------*/
/*                                                                           */
/* Program:  IMAGE.H                                                         */
/*                                                                           */
/* Purpose:  This file contains the constant and type definitions for        */
/*           the routines defined in image.c.                                */
/*                                                                           */
/* Author:   John Gauch - Version 2                                          */
/*           Zimmerman, Entenman, Fitzpatrick, Whang - Version 1             */
/*                                                                           */
/* Date:     February 23, 1987                                               */
/*                                                                           */
/*---------------------------------------------------------------------------*/

#pragma once

#ifdef SYSV
#include <sys/types.h>
#ifndef __FCNTL_HEADER__	/* Because system V does not do this */
#define __FCNTL_HEADER__	/* 	in fcntl.h		     */
#include <fcntl.h>
#endif /*!__FCNTL_HEADER__*/
#endif

#ifdef SYSV
#ifndef __FILE_HEADER__		/* Because system V does not do this */
#define __FILE_HEADER__		/* 	in file.h		     */
#include <sys/file.h>
#endif /*!__FILE_HEADER__*/
#else
  #ifdef WIN32
    //
  #else
    #include <sys/file.h>
  #endif
#endif

#ifdef WIN32
  //
#else
  #include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Boolean values */
#define TRUE		1
#define FALSE		0

/* Routine return codes */
#define VALID		1
#define INVALID		0

/* Maximum and minimum GREYTYPE values */
#define MINVAL		-32768
#define MAXVAL		32767

/* Pixel format codes */
#define GREY		0010
#define COLOR		0020
#define COLORPACKED	0040
#define USERPACKED	0200
#define BYTE		0001
#define SHORT		0002
#define LONG		0003
#define REAL		0004
#define COMPLEX		0005
#define INT             0006

/* Pixel format types */
typedef short GREYTYPE;
typedef short COLORTYPE;
typedef struct { unsigned char r,g,b,a; } CPACKEDTYPE;
typedef int USERTYPE;
typedef unsigned char BYTETYPE;
typedef short SHORTTYPE;
typedef long LONGTYPE;
typedef float REALTYPE;
typedef struct { float re, im; } COMPLEXTYPE;

/* Constants for open calls */
#define READ		(O_RDONLY)
#define UPDATE		(O_RDWR)
#define CREATE		(O_RDWR | O_CREAT | O_EXCL)

/* compression constants */
#define COMPRESS		16384
#define FORCE_COMPRESS		32768
#define FORCE_DECOMPRESS	65536

/* currently support up to 10 compression methods */
#define MAX_NUM_COMP_METHODS 10
 
/* Constants for imgetdesc calls */
#define MINMAX		0
#define HISTO		1

/* Protection modes for imcreat */
#define UOWNER		0600
#define UGROUP		0060
#define RGROUP		0040
#define UOTHER		0006
#define ROTHER		0004
#define DEFAULT		0644

/* Array length constants */
#define nADDRESS	9
#define nTITLE		81
#define nMAXMIN		2
#define nHISTOGRAM	1024
#define nDIMV		10
#define nERROR		200
#define nINFO		100

/* Old names for array length constants (from version 1) */
#define _NDIMS		nDIMV
#define _ERRLEN		nERROR
#define TITLESIZE	nTITLE
#define MAXPIX		(nHISTOGRAM-1)

/* Structure for image information (everything but pixels) */
typedef struct {
   int   Fd;			/* Computed fields */
   int   PixelSize;
   int   PixelCnt;
   int   SwapNeeded;

   int	 Compressed;		/* is pixel data compressed? */
   int	 CompressionMethod;	/* how is pixel data compressed? */
   int	 PixelsAccessed;	/* have the pixels been accessed yet? */
   int	 PixelsModified;	/* have the pixels been modified yet? */
   int	 UCPixelsFd;		/* where is the uncompressed data? */
   char	 UCPixelsFileName[256];	/* name of the uncompressed data file */

   int   Address[nADDRESS];	/* Header fields from file */
   char  Title[nTITLE];
   int   ValidMaxMin;
   int   MaxMin[nMAXMIN];
   int   ValidHistogram;
   int   Histogram[nHISTOGRAM];
   int   PixelFormat;
   int   Dimc;
   int   Dimv[nDIMV];

   int   InfoCnt;		/* Information fields from file */
   char *InfoName[nINFO];
   char *InfoData[nINFO];

   int   nImgFormat;

   } IMAGE;

/* Error string buffer */
char _imerrbuf[nERROR];

/* Initialization routines */
IMAGE *imcreat(char *new_name, int mode, int type, int dimc, int dimv[3]);
IMAGE *imopen();
IMAGE *dcmopen();
IMAGE *ifopen();
int imclose(IMAGE *);

/* Pixel access routines */
int imread();
int imwrite(IMAGE *out_image, int start, int end, float *dat);
int imgetpix();
int imputpix();

int compressImage();
int decompressImage();
int readCompressionConfigFile();
int imheaderC();

/* compression structure */
typedef struct
{
  char methodName[80];
  char compressionCommand[80];
  char decompressionCommand[80];
} compMethod;


/* Information access routines */
int imgetheader();
int imdim();
int imbounds();
int imgetdesc();
int imgettitle();
int imputtitle();
char *imgetinfo();
int imputinfo();
int imcopyinfo();
char **iminfoids();
char *imerror(); 

/*---------------------------------------------------------------------------*/
/*                                                                           */
/* File:     PDIM.H                                                          */
/*                                                                           */
/* Purpose:  This file contains declarations used by PDIM.C                  */
/*                                                                           */
/* Author:   John Gauch - Version 2                                          */
/*           Chuck Mosher - Version 1                                        */
/*                                                                           */
/* Date:     July 21, 1986                                                   */
/*                                                                           */
/*---------------------------------------------------------------------------*/

/* Global PDIM constants */
#define REC_SIZE 200

/* Types of units possible */
#define MILLIMETER 0
#define CENTIMETER 1

/* Structure for single slice description */
typedef struct {
   float Ox,Oy,Oz;
   float Ux,Uy,Uz;
   float Vx,Vy,Vz;
   float time;
   int   number;
   } SLICEREC;

/* Structure for whole PDIM description */
typedef struct {
   int version;
   int units;
   int machine;
   int slicecnt;
   SLICEREC *patient;
   SLICEREC *table;
   } PDIMREC;
 
/* Forward declarations of PDIM routines */
int pdim_read();
int pdim_write();
int pdim_free();
int pdim_append();
int pdim_window();
int pdim_scale();
int pdim_rotate();
int pdim_translate();
int pdim_map();


/* Fix up linking for the FFT library */
#ifdef IBM
#define cmplft_ cmplft
#endif


/* Fix up linking of FFT library for stardent 3000*/
#ifdef ARDENT
#define hermft_	HERMFT
#define realft_	REALFT
#define rsymft_	RSYMFT
#define sdiad_	SDIAD
#define inv21_	INV21
#define cmplft_	CMPLFT
#define srfp_	SRFP
#define diprp_	DIPRP
#define mdftkd_	MDFTKD
#define r2cftk_	R2CFTK
#define r3cftk_	R3CFTK
#define r4cftk_	R4CFTK
#define r5cftk_	R5CFTK
#define r8cftk_	R8CFTK
#define rpcftk_	RPCFTK
#endif

/******************************************************************************\
$Log: image.h,v $
 * Revision 1.3  1995/05/10  15:24:15  root
 * Revised to support compressed images
 *
\******************************************************************************/
