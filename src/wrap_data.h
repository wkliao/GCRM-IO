/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: wrap_data.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_WRAP_DATA
#define H_WRAP_DATA

#if HAVE_CONFIG_H
#include "config.h"
#endif

/*
   MODULE wrap_data
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  PURPOSE:  perform boundary updates of ghost cell along edges 
!            of subdomain blocks
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

typedef struct wrap_list_node wrap_list_node;
struct wrap_list_node {
    int tag_glbl;
    int i1;
    int j1;
    int nsd1;
    wrap_list_node *next;
};

typedef struct wrap_proc_node wrap_proc_node;
struct wrap_proc_node {
    int proc_nmbr;
    int send_total;
    int recv_total;
    int send_msgtag;
    int recv_msgtag;
    int *send_lst;
    int *i0;
    int *j0;
    int *nsd0;
    int *recv_lst;
    int *i1;
    int *j1;
    int *nsd1;
    double   **send_rk2;
    double ****send_rk4;
    double   **recv_rk2;
    double ****recv_rk4;
    wrap_list_node *send;
    wrap_list_node *recv;
    wrap_proc_node *next;
};

typedef struct wrap_node wrap_node;
struct wrap_node {
    char wrap_name[128];
    int l_wrap;
    int *l_sbdmn_pentagon_north;
    int *l_sbdmn_pentagon_south;
    int *l_sbdmn_north_pole;
    int *l_sbdmn_south_pole;
    int npe_comm;
    int rnk_comm;
    MPI_Comm comm;
    int proc_total;
    int send_message_total;
    int recv_message_total;
    MPI_Datatype flt_type;
    int sbdmn_iota;
    int *sbdmn_lst; /* pointer to an array of size [sbdmn_len] */
    int  sbdmn_len;
    int **sbdmn_perimeter_type;
    wrap_proc_node *proc;
    wrap_node *next;
};

typedef struct {
    int l_allocate_wrap_head; /* initialized to .TRUE. */
    int l_wrptmr;             /* initialized to .FALSE. */
    wrap_node *wrap_head;
    int msgtag;               /* initialized to 19 */
} MODULE_wrap_data;

#include "grid_params.h"

/* API declarations */
void init_MODULE_wrap_data(MODULE_wrap_data *wrap_data);
void initialize_wrap(MODULE_wrap_data   *wrap_data,
                     MODULE_grid_params *grid_params,
                     char            *communicator_name,
                     char            *wrap_name,
                     int              level,
                     int              subdomain_iota,
                     int             *subdomain_list,
                     int              l_report);
void finalize_MODULE_wrap_data(MODULE_wrap_data *wrap_data);

void wrap_face(MODULE_wrap_data    *wrap_data,
               const char          *wrap_name,
               double           ****face,
               const int           *dimlen);
void wrap_vrtx(MODULE_wrap_data      *wrap_data,
               const char            *wrap_name,
               double           ******vrtx,
               const int             *dimlen);
void wrap_edge(MODULE_wrap_data      *wrap_data,
               const char            *wrap_name,
               double           ******edge,
               const int             *dimlen);
void wrap_face_1lyr(MODULE_wrap_data    *wrap_data,
                    const char          *wrap_name,
                    double           ****face_1lyr,
                    const int           *dimlen3,
                    int                  idx);
void wrap_vrtx_1lyr(MODULE_wrap_data     *wrap_data,
                    const char           *wrap_name,
                    double           *****vrtx_1lyr,
                    const int            *dimlen5);
void wrap_vrtx_scalar(MODULE_wrap_data     *wrap_data,
                      const char           *wrap_name,
                      double           *****vrtx_scalar,
                      const int            *dimlen5);
void wrap_vrtx_scalar_1lyr(MODULE_wrap_data    *wrap_data,
                           const char          *wrap_name,
                           double           ****vrtx_scalar_1lyr,
                           const int           *dimlen4);
void wrap_edge_1lyr(MODULE_wrap_data      *wrap_data,
                    const char            *wrap_name,
                    double            *****edge_1lyr,
                    const int             *dimlen5);
void wrap_edge_scalar(MODULE_wrap_data     *wrap_data,
                      const char           *wrap_name,
                      double           *****edge_scalar,
                      const int            *dimlen5);
void wrap_edge_scalar_1lyr(MODULE_wrap_data    *wrap_data,
                           const char          *wrap_name,
                           double           ****edge_scalar_1lyr,
                           const int           *dimlen4);


#endif
