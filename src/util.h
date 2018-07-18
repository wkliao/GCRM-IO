/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: util.h 4597 2017-12-07 07:00:37Z wkliao $
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <mpi.h>

#define BCAST_INT(buf)     MPI_Bcast(&(buf),   1, MPI_INT,   0, MPI_COMM_WORLD);
#define BCAST_FLT(buf)     MPI_Bcast(&(buf),   1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#define BCAST_STR(buf,len) MPI_Bcast(buf,    len, MPI_CHAR,  0, MPI_COMM_WORLD);

#ifndef POWER2
#define POWER2(x) (1 << (x))
#endif
#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#define OPEN_FILE(fname) {                                            \
    fp = fopen(fname, "r");                                           \
    if (fp == NULL) {                                                 \
        printf("Error (file=%s line=%d): file open %s (%s)\n",        \
               __FILE__, __LINE__, fname, strerror(errno));           \
        MPI_Abort(MPI_COMM_WORLD, -1);                                \
    }                                                                 \
}

#define GET_LINE_TOKEN_PAIR {                                     \
    /* line contains a '\n' at the end */                         \
                                                                  \
    key = strtok(line, "= \t\n");                                 \
    if (key == NULL || key[0] == '#' || key[0] == '!') continue;  \
                                                                  \
    value = strtok(NULL, " !\t\n");                               \
}


#define ABORT {printf("MPI_Abort: at file=%s line=%d func=%s\n", __FILE__,__LINE__,__func__); fflush(stdout); MPI_Abort(MPI_COMM_WORLD, -1);}


/* Check the status return value from a io call.
   message string a string to associate with error
   nodie logical override the default abort behavior which is useful in cases
   such as opening a file that may not exist.  Default is false which means
   an error will NOT lead to an abort.
*/
#define NC_CHECK(status, message, nodie)                                      \
    if (status != NC_NOERR) {                                                 \
        printf("Error: %s %d %s (at file=%s line=%d func=%s)\n", message,     \
               status, ncmpi_strerror(status), __FILE__,__LINE__,__func__);   \
        if (nodie) ABORT                                                      \
    }

#define tmalloc(a)    malloc_fn(a,__LINE__,__func__,__FILE__)
#define tcalloc(a,b)  calloc_fn(a,b,__LINE__,__func__,__FILE__)
#define tfree(a)      free_fn(a,__LINE__,__func__,__FILE__)

void* malloc_fn(size_t size, const int lineno, const char *func, const char *filename);

void* calloc_fn(size_t count, size_t size, const int lineno, const char *func, const char *filename);

void free_fn(void *buf, const int lineno, const char *func, const char *filename);

long long get_max_mem_alloc(void);
long long get_mem_alloc(void);
void reset_mem_alloc(void);
void check_mem_alloc(void);

int    * calloc_1D_int(size_t x);
float  * calloc_1D_flt(size_t x);
double * calloc_1D_dbl(size_t x);
int    ** calloc_1D_intptr(size_t x);

void   free_1D_int(int    *buf);
void   free_1D_flt(float  *buf);
void   free_1D_dbl(double *buf);
void   free_1D_intptr(int **buf);

int    ** calloc_2D_int(size_t y, size_t x);
float  ** calloc_2D_flt(size_t y, size_t x);
double ** calloc_2D_dbl(size_t y, size_t x);
int    *** calloc_2D_intptr(size_t y, size_t x);
double *** calloc_2D_dblptr(size_t y, size_t x);

void   free_2D_int(int    **buf);
void   free_2D_flt(float  **buf);
void   free_2D_dbl(double **buf);
void   free_2D_intptr(int ***buf);
void   free_2D_dblptr(double ***buf);

int    *** calloc_3D_int(size_t z, size_t y, size_t x);
float  *** calloc_3D_flt(size_t z, size_t y, size_t x);
double *** calloc_3D_dbl(size_t z, size_t y, size_t x);

void   free_3D_int(int    ***buf);
void   free_3D_flt(float  ***buf);
void   free_3D_dbl(double ***buf);

int    **** calloc_4D_int(size_t z, size_t y, size_t x, size_t w);
float  **** calloc_4D_flt(size_t z, size_t y, size_t x, size_t w);
double **** calloc_4D_dbl(size_t z, size_t y, size_t x, size_t w);
double ***** calloc_4D_dblptr(size_t z, size_t y, size_t x, size_t w);

void   free_4D_int(int    ****buf);
void   free_4D_flt(float  ****buf);
void   free_4D_dbl(double ****buf);
void   free_4D_dblptr(double *****buf);

int    ***** calloc_5D_int(size_t z, size_t y, size_t x, size_t w, size_t v);
float  ***** calloc_5D_flt(size_t z, size_t y, size_t x, size_t w, size_t v);
double ***** calloc_5D_dbl(size_t z, size_t y, size_t x, size_t w, size_t v);

void   free_5D_int(int    *****buf);
void   free_5D_flt(float  *****buf);
void   free_5D_dbl(double *****buf);

int    ****** calloc_6D_int(size_t z, size_t y, size_t x, size_t w, size_t v, size_t u);
float  ****** calloc_6D_flt(size_t z, size_t y, size_t x, size_t w, size_t v, size_t u);
double ****** calloc_6D_dbl(size_t z, size_t y, size_t x, size_t w, size_t v, size_t u);

void   free_6D_int(int    ******buf);
void   free_6D_flt(float  ******buf);
void   free_6D_dbl(double ******buf);

void random_1D_int(int         *buf, int z);
void random_1D_flt(float       *buf, int z);
void random_1D_dbl(double      *buf, int z);

void random_2D_int(int        **buf, int z, int y);
void random_2D_flt(float      **buf, int z, int y);
void random_2D_dbl(double     **buf, int z, int y);

void random_3D_int(int       ***buf, int z, int y, int x);
void random_3D_flt(float     ***buf, int z, int y, int x);
void random_3D_dbl(double    ***buf, int z, int y, int x);

void random_4D_int(int      ****buf, int z, int y, int x, int w);
void random_4D_flt(float    ****buf, int z, int y, int x, int w);
void random_4D_dbl(double   ****buf, int z, int y, int x, int w);

void random_5D_int(int     *****buf, int z, int y, int x, int w, int v);
void random_5D_flt(float   *****buf, int z, int y, int x, int w, int v);
void random_5D_dbl(double  *****buf, int z, int y, int x, int w, int v);

void random_6D_int(int    ******buf, int z, int y, int x, int w, int v, int u);
void random_6D_flt(float  ******buf, int z, int y, int x, int w, int v, int u);
void random_6D_dbl(double ******buf, int z, int y, int x, int w, int v, int u);

void cshift_1D(void *buf, int x, int elm_size, int shift);
void cshift_2D(void *buf, int y, int x, int elm_size, int shift);
void cshift_down_2D(void *buf, int n, int dim_size);
void cshift_up_2D(void *buf, int n, int dim_size);

void print_mpi_info(MPI_Info *info_used);

#endif
