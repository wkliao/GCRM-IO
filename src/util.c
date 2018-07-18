/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: util.c 4597 2017-12-07 07:00:37Z wkliao $
 */

#include <stdio.h>
#include <stdlib.h>
#include <search.h> /* tsearch() and tdelete() */
#include <string.h> /* memcpy() memset() */
#include <mpi.h>
#include "util.h"

/* to enable malloc tracking */
#define TRACK_MALLOC

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

static long long mem_alloc, max_mem_alloc;

#ifdef USE_LINEAR_TABLE_LOOKUP
#define USE_LINEAR_TABLE_LOOKUP
typedef struct {
    void   *buf;
    size_t  size;
} mem_entry;

static mem_entry *mem_table;
static int alloc_len;
static int cur_len;
static int tail_index;

#else
typedef struct {
    void   *self;
    void   *buf;
    size_t  size;
    int     lineno;
    char   *func;
    char   *filename;
} mem_entry;

static void *mem_root;

static
int cmp(const void *a, const void *b) {
    mem_entry *fa = (mem_entry*)a;
    mem_entry *fb = (mem_entry*)b;

    if (fa->buf > fb->buf) return  1;
    if (fa->buf < fb->buf) return -1;
    return 0;
}

static
void walker(const void *node, const VISIT which, const int depth) {
    mem_entry *f;
    f = *(mem_entry **)node;
    printf("Warning: malloc yet to be freed (buf=%p size=%zd filename=%s func=%s line=%d)\n", f->buf, f->size, f->filename, f->func, f->lineno);

}
#endif

/*----< reset_mem_alloc() >---------------------------------------------------*/
void reset_mem_alloc(void)
{
        mem_alloc = 0;
    max_mem_alloc = 0;

#ifdef USE_LINEAR_TABLE_LOOKUP
    tail_index = 0;
    cur_len = 0;
    alloc_len = 1024;
    mem_table = (mem_entry*) calloc(alloc_len, sizeof(mem_entry));
#else
    mem_root = NULL;
#endif
}

/*----< add_mem_entry() >-----------------------------------------------------*/
/* add a new malloc entry to the table */
static
void add_mem_entry(void   *buf,
                   size_t  size,
                   const int   lineno,
                   const char *func,
                   const char *filename)
{
#ifndef TRACK_MALLOC
    return;
#endif

#ifdef USE_LINEAR_TABLE_LOOKUP
    int indx;

    if (cur_len == alloc_len) { /* expand the table size */
        alloc_len += 1024;
        mem_table = (mem_entry*) realloc(mem_table, alloc_len * sizeof(mem_entry));
        memset(mem_table+cur_len, 0, 1024*sizeof(mem_entry));
    }

    if (cur_len == tail_index) { /* no empty slot before tail */
        indx = cur_len;
        tail_index++;
    }
    else {  /* find the empty slot and fill in */
        for (indx=tail_index-1; indx>=0; indx--) {
            if (mem_table[indx].buf == NULL)
                break;
        }
        if (indx < 0) {
            printf("Error: add_mem_entry()\n");
            return;
        }
    }
    mem_table[indx].buf  = buf;
    mem_table[indx].size = size;
    cur_len++;

#else
    /* use C tsearch utility */
    mem_entry *node = (mem_entry*) malloc(sizeof(mem_entry));
    node->self     = node;
    node->buf      = buf;
    node->size     = size;
    node->lineno   = lineno;
    node->func     = (char*)malloc(strlen(func)+1);
    node->filename = (char*)malloc(strlen(filename)+1);
    strcpy(node->func, func);
    node->func[strlen(func)] = '\0';
    strcpy(node->filename, filename);
    node->filename[strlen(filename)] = '\0';

    void *ret = tsearch(node, &mem_root, cmp);
    if (ret == NULL) {
        fprintf(stderr, "Error at line %d file %s: tsearch()\n",
                __LINE__,__FILE__);
        return;
    }
#endif
    mem_alloc += size;
    max_mem_alloc = MAX(max_mem_alloc, mem_alloc);
}

/*----< del_mem_entry() >-----------------------------------------------------*/
/* delete a malloc entry from the table */
static
void del_mem_entry(void *buf)
{
#ifndef TRACK_MALLOC
    return;
#endif

#ifdef USE_LINEAR_TABLE_LOOKUP
    int indx;

    for (indx=tail_index-1; indx>=0; indx--) {
        if (mem_table[indx].buf == buf)
            break;
    }
    if (indx < 0) {
        printf("Error: del_mem_entry()\n");
        return;
    }

    mem_alloc -= mem_table[indx].size;
    mem_table[indx].buf  = NULL;
    mem_table[indx].size = 0;
    if (indx == tail_index-1) {
        while (mem_table[indx].buf == NULL && indx >=0)
            indx--;
        tail_index = indx + 1;
    }
    cur_len--;
#else
    /* use C tsearch utility */
    if (mem_root != NULL) {
        mem_entry node;
        node.buf  = buf;
        void *ret = tfind(&node, &mem_root, cmp);
        mem_entry **found = (mem_entry**) ret;
        if (ret == NULL) {
            fprintf(stderr, "Error at line %d file %s: tfind() buf=%p\n",
                    __LINE__,__FILE__,buf);
            return;
        }
        /* free space for func and filename */
        free((*found)->func);
        free((*found)->filename);

        /* subtract the space amount to be freed */
        mem_alloc -= (*found)->size;
        void *tmp = (*found)->self;
        ret = tdelete(&node, &mem_root, cmp);
        if (ret == NULL) {
            fprintf(stderr, "Error at line %d file %s: tdelete() buf=%p\n",
                    __LINE__,__FILE__,buf);
            return;
        }
        free(tmp);
    }
    else
        fprintf(stderr, "Error at line %d file %s: mem_root is NULL\n",
                __LINE__,__FILE__);
#endif
}

/*----< get_mem_alloc() >-----------------------------------------------------*/
/* get the current aggregate size allocated by malloc */
long long get_mem_alloc(void)
{
    return mem_alloc;
}

/*----< get_max_mem_alloc() >-------------------------------------------------*/
/* get the max watermark ever researched by malloc */
long long get_max_mem_alloc(void)
{
    return max_mem_alloc;
}

/*----< check_mem_alloc() >---------------------------------------------------*/
/* check if there is any malloc residue */
void check_mem_alloc(void)
{
#ifdef USE_LINEAR_TABLE_LOOKUP
    /* get the current number of malloc calls made, but have not freed yet */
    if (cur_len > 0)
        printf("Warning: memory allocated yet to be freed num_cur = %d\n", cur_len);

    /* get the tail of the malloc register table */
    if (tail_index > 0)
        printf("Warning: memory allocated yet to be freed tail = %d\n", tail_index);
#else
    /* check if malloc tree is empty */
    if (mem_root != NULL) {
        printf("Warning: there are malloc residues\n");
        twalk(mem_root, walker);
    }
#endif
}

#define ERR_CHECK(buf, len) if (len > 0 && buf == NULL) { printf("Error: malloc failed len=%zd (file=%s line=%d)\n", len, __FILE__, __LINE__); ABORT }

#define tmalloc(a)    malloc_fn(a,__LINE__,__func__,__FILE__)
#define tcalloc(a,b)  calloc_fn(a,b,__LINE__,__func__,__FILE__)
#define tfree(a)      free_fn(a,__LINE__,__func__,__FILE__)

/*----< malloc_fn() >-----------------------------------------------------------*/
void* malloc_fn(size_t      size,
                const int   lineno,
                const char *func,
                const char *filename)
{
    void *buf = malloc(size);
    ERR_CHECK(buf, size)
    add_mem_entry(buf, size, lineno, func, filename);
    return buf;
}

/*----< calloc_fn() >-----------------------------------------------------------*/
void* calloc_fn(size_t count,
                size_t size,
                const int   lineno,
                const char *func,
                const char *filename)
{
    void *buf = calloc(count, size);
    ERR_CHECK(buf, count*size)
    add_mem_entry(buf, count * size, lineno, func, filename);
    return buf;
}

/*----< free_fn() >-------------------------------------------------------------*/
void free_fn(void *buf,
             const int   lineno,
             const char *func,
             const char *filename)
{
    del_mem_entry(buf);
    free(buf);
}

#define CALLOC_1D(fname, type)                      \
type * calloc_1D_##fname(size_t x)                  \
{                                                   \
    size_t req_len = x * sizeof(type);              \
    type *buf = (type*) calloc(x, sizeof(type));    \
    ERR_CHECK(buf, x)                               \
    add_mem_entry(buf, req_len, __LINE__, __func__, __FILE__); \
    return buf;                                     \
}

CALLOC_1D(dbl, double)
CALLOC_1D(flt, float)
CALLOC_1D(int, int)

#define FREE_1D(fname, type)                      \
void free_1D_##fname(type *buf)                   \
{                                                 \
    del_mem_entry(buf);                           \
    free(buf);                                    \
}

int ** calloc_1D_intptr(size_t x)
{
    size_t req_len = x * sizeof(int*);
    int **buf = (int**) calloc(x, sizeof(int*));
    ERR_CHECK(buf, x)
    add_mem_entry(buf, req_len, __LINE__, __func__, __FILE__);
    return buf;
}

void free_1D_intptr(int **buf)
{
    del_mem_entry(buf);
    free(buf);
}

FREE_1D(dbl, double)
FREE_1D(flt, float)
FREE_1D(int, int)

#define CALLOC_2D(fname, type)                          \
type ** calloc_2D_##fname(size_t y, size_t x)           \
{                                                       \
    int _j;                                             \
    size_t req_len = y*  sizeof(type*)                  \
                   + y*x*sizeof(type);                  \
    type **buf  = (type**) malloc(y   * sizeof(type*)); \
    ERR_CHECK(buf, y)                                   \
    type  *bufx = (type*)  calloc(y*x,  sizeof(type));  \
    ERR_CHECK(bufx, y*x)                                \
    for (_j=0; _j<y; _j++, bufx+=x)                     \
        buf[_j] = bufx;                                 \
    add_mem_entry(buf, req_len, __LINE__, __func__, __FILE__); \
    return buf;                                         \
}

CALLOC_2D(dbl, double)
CALLOC_2D(flt, float)
CALLOC_2D(int, int)
CALLOC_2D(dblptr, double*)
CALLOC_2D(intptr, int*)

#define FREE_2D(fname, type)                       \
void free_2D_##fname(type **buf)                   \
{                                                  \
    del_mem_entry(buf);                            \
    free(buf[0]);                                  \
    free(buf);                                     \
}

FREE_2D(dbl, double)
FREE_2D(flt, float)
FREE_2D(int, int)
FREE_2D(dblptr, double*)
FREE_2D(intptr, int*)

#define CALLOC_3D(fname, type)                               \
type *** calloc_3D_##fname(size_t z, size_t y, size_t x)     \
{                                                            \
    int _j, _k;                                              \
    size_t req_len = z*    sizeof(type**)                    \
                   + z*y*  sizeof(type*)                     \
                   + z*y*x*sizeof(type);                     \
    type ***buf  = (type***) malloc(z     * sizeof(type**)); \
    ERR_CHECK(buf, z)                                        \
    type  **bufy = (type**)  malloc(z*y   * sizeof(type*));  \
    ERR_CHECK(bufy, z*y)                                     \
    type   *bufx = (type*)   calloc(z*y*x,  sizeof(type));   \
    ERR_CHECK(bufx, z*y*x)                                   \
    for (_k=0; _k<z; _k++, bufy+=y) {                        \
        buf[_k] = bufy;                                      \
        for (_j=0; _j<y; _j++, bufx+=x)                      \
            buf[_k][_j] = bufx;                              \
    }                                                        \
    add_mem_entry(buf, req_len, __LINE__, __func__, __FILE__); \
    return buf;                                              \
}


CALLOC_3D(dbl, double)
CALLOC_3D(flt, float)
CALLOC_3D(int, int)

#define FREE_3D(fname, type)                       \
void free_3D_##fname(type ***buf)                  \
{                                                  \
    del_mem_entry(buf);                            \
    free(buf[0][0]);                               \
    free(buf[0]);                                  \
    free(buf);                                     \
}

void free_3D_dbl(double ***buf)
{
    del_mem_entry(buf);
    free(buf[0][0]);  
    free(buf[0]);    
    free(buf);      
}

// FREE_3D(dbl, double)
FREE_3D(flt, float)
FREE_3D(int, int)

#define CALLOC_4D(fname, type)                                      \
type **** calloc_4D_##fname(size_t z, size_t y, size_t x, size_t w) \
{                                                                   \
    int _i, _j, _k;                                                 \
    size_t req_len = z*      sizeof(type***)                        \
                   + z*y*    sizeof(type**)                         \
                   + z*y*x*  sizeof(type*)                          \
                   + z*y*x*w*sizeof(type);                          \
    type ****buf  = (type****) malloc(z       * sizeof(type***));   \
    ERR_CHECK(buf, z)                                               \
    type  ***bufz = (type***)  malloc(z*y     * sizeof(type**));    \
    ERR_CHECK(bufz, z*y)                                            \
    type   **bufy = (type**)   malloc(z*y*x   * sizeof(type*));     \
    ERR_CHECK(bufy, z*y*x)                                          \
    type    *bufx = (type*)    calloc(z*y*x*w,  sizeof(type));      \
    ERR_CHECK(bufx, z*y*x*w)                                        \
    for (_k=0; _k<z; _k++, bufz+=y) {                               \
        buf[_k] = bufz;                                             \
        for (_j=0; _j<y; _j++, bufy+=x) {                           \
            buf[_k][_j] = bufy;                                     \
            for (_i=0; _i<x; _i++, bufx+=w)                         \
                buf[_k][_j][_i] = bufx;                             \
        }                                                           \
    }                                                               \
    add_mem_entry(buf, req_len, __LINE__, __func__, __FILE__);      \
    return buf;                                                     \
}

CALLOC_4D(dbl, double)
CALLOC_4D(flt, float)
CALLOC_4D(int, int)
CALLOC_4D(dblptr, double*)

#define FREE_4D(fname, type)                       \
void free_4D_##fname(type ****buf)                 \
{                                                  \
    del_mem_entry(buf);                            \
    free(buf[0][0][0]);                            \
    free(buf[0][0]);                               \
    free(buf[0]);                                  \
    free(buf);                                     \
}

FREE_4D(dbl, double)
FREE_4D(flt, float)
FREE_4D(int, int)
FREE_4D(dblptr, double*)

#define CALLOC_5D(fname, type)                                                 \
type ***** calloc_5D_##fname(size_t z, size_t y, size_t x, size_t w, size_t v) \
{                                                                              \
    int _i, _j, _k, _l;                                                        \
    size_t req_len = z*        sizeof(type****)                                \
                   + z*y*      sizeof(type***)                                 \
                   + z*y*x*    sizeof(type**)                                  \
                   + z*y*x*w*  sizeof(type*)                                   \
                   + z*y*x*w*v*sizeof(type);                                   \
    type *****buf  = (type*****) malloc(z         * sizeof(type****));         \
    ERR_CHECK(buf, z)                                                          \
    type  ****bufz = (type****)  malloc(z*y       * sizeof(type***));          \
    ERR_CHECK(bufz, z*y)                                                       \
    type   ***bufy = (type***)   malloc(z*y*x     * sizeof(type**));           \
    ERR_CHECK(bufy, z*y*x)                                                     \
    type    **bufx = (type**)    malloc(z*y*x*w   * sizeof(type*));            \
    ERR_CHECK(bufx, z*y*x*w)                                                   \
    type     *bufw = (type*)     calloc(z*y*x*w*v,  sizeof(type));             \
    ERR_CHECK(bufw, z*y*x*w*v)                                                 \
    for (_l=0; _l<z; _l++, bufz+=y) {                                          \
        buf[_l] = bufz;                                                        \
        for (_k=0; _k<y; _k++, bufy+=x) {                                      \
            buf[_l][_k] = bufy;                                                \
            for (_j=0; _j<x; _j++, bufx+=w) {                                  \
                buf[_l][_k][_j] = bufx;                                        \
                for (_i=0; _i<w; _i++, bufw+=v)                                \
                    buf[_l][_k][_j][_i] = bufw;                                \
            }                                                                  \
        }                                                                      \
    }                                                                          \
    add_mem_entry(buf, req_len, __LINE__, __func__, __FILE__);                 \
    return buf;                                                                \
}

CALLOC_5D(dbl, double)
CALLOC_5D(flt, float)
CALLOC_5D(int, int)

#define FREE_5D(fname, type)                       \
void free_5D_##fname(type *****buf)                \
{                                                  \
    del_mem_entry(buf);                            \
    free(buf[0][0][0][0]);                         \
    free(buf[0][0][0]);                            \
    free(buf[0][0]);                               \
    free(buf[0]);                                  \
    free(buf);                                     \
}

FREE_5D(dbl, double)
FREE_5D(flt, float)
FREE_5D(int, int)

#define CALLOC_6D(fname, type)                                                            \
type ****** calloc_6D_##fname(size_t z, size_t y, size_t x, size_t w, size_t v, size_t u) \
{                                                                                         \
    int _i, _j, _k, _l, _m;                                                               \
    size_t req_len = z*          sizeof(type*****)                                        \
                   + z*y*        sizeof(type****)                                         \
                   + z*y*x*      sizeof(type***)                                          \
                   + z*y*x*w*    sizeof(type**)                                           \
                   + z*y*x*w*v*  sizeof(type*)                                            \
                   + z*y*x*w*v*u*sizeof(type);                                            \
    type ******buf  = (type******) malloc(z           * sizeof(type*****));               \
    ERR_CHECK(buf, z)                                                                     \
    type  *****bufz = (type*****)  malloc(z*y         * sizeof(type****));                \
    ERR_CHECK(bufz, z*y)                                                                  \
    type   ****bufy = (type****)   malloc(z*y*x       * sizeof(type***));                 \
    ERR_CHECK(bufy, z*y*x)                                                                \
    type    ***bufx = (type***)    malloc(z*y*x*w     * sizeof(type**));                  \
    ERR_CHECK(bufx, z*y*x*w)                                                              \
    type     **bufw = (type**)     malloc(z*y*x*w*v   * sizeof(type*));                   \
    ERR_CHECK(bufw, z*y*x*w*v)                                                            \
    type      *bufv = (type*)      calloc(z*y*x*w*v*u,  sizeof(type));                    \
    ERR_CHECK(bufv, z*y*x*w*v*u)                                                          \
    for (_m=0; _m<z; _m++, bufz+=y) {                                                     \
        buf[_m] = bufz;                                                                   \
        for (_l=0; _l<y; _l++, bufy+=x) {                                                 \
            buf[_m][_l] = bufy;                                                           \
            for (_k=0; _k<x; _k++, bufx+=w) {                                             \
                buf[_m][_l][_k] = bufx;                                                   \
                for (_j=0; _j<w; _j++, bufw+=v) {                                         \
                    buf[_m][_l][_k][_j] = bufw;                                           \
                    for (_i=0; _i<v; _i++, bufv+=u)                                       \
                        buf[_m][_l][_k][_j][_i] = bufv;                                   \
                }                                                                         \
            }                                                                             \
        }                                                                                 \
    }                                                                                     \
    add_mem_entry(buf, req_len, __LINE__, __func__, __FILE__);                            \
    return buf;                                                                           \
}

CALLOC_6D(dbl, double)
CALLOC_6D(flt, float)
CALLOC_6D(int, int)

#define FREE_6D(fname, type)                       \
void free_6D_##fname(type ******buf)               \
{                                                  \
    del_mem_entry(buf);                            \
    free(buf[0][0][0][0][0]);                      \
    free(buf[0][0][0][0]);                         \
    free(buf[0][0][0]);                            \
    free(buf[0][0]);                               \
    free(buf[0]);                                  \
    free(buf);                                     \
}

FREE_6D(dbl, double)
FREE_6D(flt, float)
FREE_6D(int, int)


#define RANDOM_1D(fname, type)                                                 \
void random_1D_##fname(type *buf, int x)                                       \
{                                                                              \
    int i;                                                                     \
    for (i=0; i<x; i++)                                                        \
        buf[i] = random();                                                     \
}

RANDOM_1D(dbl, double)
RANDOM_1D(flt, float)
RANDOM_1D(int, int)

#define RANDOM_2D(fname, type)                                                 \
void random_2D_##fname(type **buf, int y, int x)                               \
{                                                                              \
    int i, j;                                                                  \
    for (i=0; i<y; i++)                                                        \
        for (j=0; j<x; j++)                                                    \
            buf[i][j] = random();                                              \
}

RANDOM_2D(dbl, double)
RANDOM_2D(flt, float)
RANDOM_2D(int, int)

#define RANDOM_3D(fname, type)                                                 \
void random_3D_##fname(type ***buf, int z, int y, int x)                       \
{                                                                              \
    int i, j, k;                                                               \
    for (i=0; i<z; i++)                                                        \
        for (j=0; j<y; j++)                                                    \
            for (k=0; k<x; k++)                                                \
                buf[i][j][k] = random();                                       \
}

RANDOM_3D(dbl, double)
RANDOM_3D(flt, float)
RANDOM_3D(int, int)

#define RANDOM_4D(fname, type)                                                 \
void random_4D_##fname(type ****buf, int z, int y, int x, int w)               \
{                                                                              \
    int i, j, k, l;                                                            \
    for (i=0; i<z; i++)                                                        \
        for (j=0; j<y; j++)                                                    \
            for (k=0; k<x; k++)                                                \
                for (l=0; l<w; l++)                                            \
                    buf[i][j][k][l] = random();                                \
}

RANDOM_4D(dbl, double)
RANDOM_4D(flt, float)
RANDOM_4D(int, int)

#define RANDOM_5D(fname, type)                                                 \
void random_5D_##fname(type *****buf, int z, int y, int x, int w, int v)       \
{                                                                              \
    int i, j, k, l, m;                                                         \
    for (i=0; i<z; i++)                                                        \
        for (j=0; j<y; j++)                                                    \
            for (k=0; k<x; k++)                                                \
                for (l=0; l<w; l++)                                            \
                    for (m=0; m<v; m++)                                        \
                        buf[i][j][k][l][m] = random();                         \
}

RANDOM_5D(dbl, double)
RANDOM_5D(flt, float)
RANDOM_5D(int, int)

#define RANDOM_6D(fname, type)                                                 \
void random_6D_##fname(type ******buf, int z, int y, int x, int w, int v, int u)\
{                                                                              \
    int i, j, k, l, m, n;                                                      \
    for (i=0; i<z; i++)                                                        \
        for (j=0; j<y; j++)                                                    \
            for (k=0; k<x; k++)                                                \
                for (l=0; l<w; l++)                                            \
                    for (m=0; m<v; m++)                                        \
                        for (n=0; n<u; n++)                                    \
                            buf[i][j][k][l][m][n] = random();                  \
}

RANDOM_6D(dbl, double)
RANDOM_6D(flt, float)
RANDOM_6D(int, int)


/*----< cshift_1D() >---------------------------------------------------------*/
void cshift_1D(void *buf,      /* [x] */
               int   x,
               int   elm_size,
               int   shift)
/* Before: V = {1, 2, 3, 4, 5, 6};
   cshift_1D(V, 6, 2, sizeof(int));
   After:  V = {3, 4, 5, 6, 1, 2};

   Before: V = {1, 2, 3, 4, 5, 6};
   cshift_1D(V, 6, -2, sizeof(int));
   After:  V = {5, 6, 1, 2, 3, 4};
*/
{
    cshift_2D(buf, x, 1, elm_size, shift);

/* below is the older version for buf of int type
    int remain, *tmp;
    if (shift == 0) return;

    if (shift > 0) {
        tmp = (int*) malloc(shift * sizeof(int));
        remain = dim_size - shift;
        memcpy(tmp, buf, shift * sizeof(int));
        memmove(buf, buf+shift, remain * sizeof(int));
        memcpy(buf+remain, tmp, shift * sizeof(int));
    }
    else {
        shift = -shift;
        tmp = (int*) malloc(shift * sizeof(int));
        remain = dim_size - shift;
        memcpy(tmp, buf+remain, shift * sizeof(int));
        memmove(buf+shift, buf, remain * sizeof(int));
        memcpy(buf, tmp, shift * sizeof(int));
    }
    free(tmp);
*/
}

/*----< cshift_2D() >---------------------------------------------------------*/
void cshift_2D(void *buf,       /* &[y][x] */
               int   y,
               int   x,
               int   elm_size,
               int   shift)
{
    int len, disp, abs;
    char *cbuf, *tmp;

    if (shift == 0) return;

    abs  = (shift > 0) ? shift : -shift;
    disp = (y-abs)*x*elm_size;
    len  = abs*x*elm_size;
    tmp  = (char*) malloc(len);
    cbuf = (char*) buf;
    if (shift > 0) {
        memcpy(tmp, cbuf, len);
        memmove(cbuf, cbuf+len, disp);
        memcpy(cbuf+disp, tmp, len);
    } else {
        memcpy(tmp, cbuf+disp, len);
        memmove(cbuf+len, cbuf, disp);
        memcpy(cbuf, tmp, len);
    }
    free(tmp);
}

/*----< cshift_up_2D() >------------------------------------------------------*/
void cshift_up_2D(void *buf,    /* &[y][x] */
                  int y,
                  int dim_size) /* x * element_size */
/* CSHIFT (buf(:,:), SHIFT=1, DIM=2) */
{
    char *cbuf = (char*) buf;
    char *tmp = (char*) malloc(dim_size);
    memcpy(tmp, cbuf, dim_size);
    memmove(cbuf, cbuf+dim_size, (y-1)*dim_size);
    memcpy(cbuf+(y-1)*dim_size, tmp, dim_size);
    free(tmp);
}

/*----< cshift_down_2D() >----------------------------------------------------*/
void cshift_down_2D(void *buf,      /* &[y][x] */
                    int   y,
                    int   dim_size) /* x * element_size */
/* CSHIFT (buf(:,:), SHIFT=-1, DIM=2) */
{
    char *cbuf = (char*) buf;
    char *tmp = (char*) malloc(dim_size);
    memcpy(tmp, cbuf+(y-1)*dim_size, dim_size);
    memmove(cbuf+dim_size, cbuf, (y-1)*dim_size);
    memcpy(cbuf, tmp, dim_size);
    free(tmp);
}

/*----< print_mpi_info() >----------------------------------------------------*/
void print_mpi_info(MPI_Info *info_used)
{
    int  i, nkeys;

    MPI_Info_get_nkeys(*info_used, &nkeys);
    printf("---- MPI file info used ----\n");
    for (i=0; i<nkeys; i++) {
        char key[MPI_MAX_INFO_KEY], value[MPI_MAX_INFO_VAL];
        int  valuelen, flag;

        MPI_Info_get_nthkey(*info_used, i, key);
        MPI_Info_get_valuelen(*info_used, key, &valuelen, &flag);
        MPI_Info_get(*info_used, key, valuelen+1, value, &flag);
        printf("MPI File Info: [%2d] key = %25s, value = %s\n",i,key,value);
    }
}

#ifdef TEST_ALONE
#define TEST_ALONE
int main(int argc, char* argv[]) {

    reset_mem_alloc();    
    int     *ibuf = calloc_1D_int(10);
    double **dbuf = calloc_2D_dbl(10, 10);
    free_2D_dbl(dbuf);
    free_1D_int(ibuf);

    printf("sizeof(double*) is %ld\n", sizeof(double*));
    printf("sizeof(int*)    is %ld\n", sizeof(int*));
    printf("max mem_alloc is %lld\n", get_max_mem_alloc());
    printf("    mem_alloc is %lld\n", get_mem_alloc());
    check_mem_alloc();
    return 0;
}
#endif
