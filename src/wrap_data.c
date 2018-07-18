/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: wrap_data.c 4608 2017-12-07 07:22:26Z wkliao $
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* memset() */
#include <mpi.h>

#include "wrap_data.h"
#include "gcrm.h"

static void set_sbdmn_perimeter_type(wrap_node *wrp,
                              grid_node *grid,
                              int        level_max,
                              int        level,
                              int        iota,
                              int       *subdomain_lst,
                              int        sbdmn_len);
static void set_recv_lst(wrap_node *wrp,
                  grid_node *grid,
                  int        level_max,
                  int        level,
                  int        iota,
                  int       *subdomain_lst,
                  int        sbdmn_len);
static wrap_proc_node* get_wrap_proc(wrap_node *wrp, int proc_nmbr);
static wrap_list_node* get_wrap_list(char *option, wrap_proc_node *proc);
static void wrap_4d(wrap_node *wrp, double ****x1, const int *dimlen);
static void wrap_6d(wrap_node *wrp, double ******x2, const int *dimlen);
static void adjust_face(wrap_node *wrp, double ****face, const int *dimlen);
static void adjust_edge(wrap_node *wrp, double ******edge, const int *dimlen);
static void adjust_vrtx(wrap_node *wrp, double ******vrtx, const int *dimlen);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  PURPOSE:  perform boundary updates of ghost cell along edges 
!            of subdomain blocks
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/*----< init_MODULE_wrap_data() >---------------------------------------------*/
void init_MODULE_wrap_data(MODULE_wrap_data *wrap_data)
{
    wrap_data->l_allocate_wrap_head = 1;
    wrap_data->l_wrptmr             = 0;
    wrap_data->msgtag               = 19;
    wrap_data->wrap_head            = NULL;
}


/*----< finalize_MODULE_wrap_data() >-----------------------------------------*/
void finalize_MODULE_wrap_data(MODULE_wrap_data *wrap_data)
{
    while (wrap_data->wrap_head != NULL) {
        wrap_node *wrp = wrap_data->wrap_head;
        wrap_data->wrap_head = wrp->next;

        free_1D_int(wrp->l_sbdmn_south_pole);
        free_1D_int(wrp->l_sbdmn_north_pole);
        free_1D_int(wrp->l_sbdmn_pentagon_south);
        free_1D_int(wrp->l_sbdmn_pentagon_north);
        free_2D_int(wrp->sbdmn_perimeter_type);

        while (wrp->proc != NULL) {
            wrap_proc_node *proc = wrp->proc;
            wrp->proc = proc->next;

            if (proc->recv_total > 0) tfree(proc->recv_lst);
            if (wrp->proc_total  > 0) tfree(proc->send_lst);
            if (proc->send_total > 0) tfree(proc->i0);

            while (proc->send != NULL) {
                wrap_list_node *instr = proc->send;
                proc->send = instr->next;
                tfree(instr);
            }
            while (proc->recv != NULL) {
                wrap_list_node *instr = proc->recv;
                proc->recv = instr->next;
                tfree(instr);
            }
            tfree(proc);
        }
        tfree(wrp);
    }
}


/*----< initialize_wrap() >---------------------------------------------------*/
void initialize_wrap(MODULE_wrap_data   *wrap_data,
                     MODULE_grid_params *grid_params,
                     char               *communicator_name,
                     char               *wrap_name,
                     int                 level,
                     int                 subdomain_iota,
                     int                *subdomain_lst,
                     int                 l_report)
{
    int sbdmn_len, n, m, proc_total, ix[3];
    wrap_node *wrp;
    wrap_proc_node *proc;
    wrap_list_node *instr;

/*-----------------------------------------------------------------------
!  the first time initialize_wrap is called allocate the head node
!-----------------------------------------------------------------------*/
    if (wrap_data->l_allocate_wrap_head) {
        wrap_data->l_allocate_wrap_head = 0;
        wrap_data->wrap_head = (wrap_node*) tmalloc(sizeof(wrap_node));
        wrp = wrap_data->wrap_head;
        strcpy(wrp->wrap_name, wrap_name);
        wrp->l_wrap = 1;
        wrp->proc_total = 0;
        wrp->send_message_total = 0;
        wrp->recv_message_total = 0;
        wrp->proc = NULL;
        wrp->next = NULL;                     
    }
    else {
/*-----------------------------------------------------------------------
!  look through the list to find duplicate wrap_name
!-----------------------------------------------------------------------*/
        wrp = wrap_data->wrap_head;
        while (wrp != NULL) {
            if (! strcmp(wrp->wrap_name, wrap_name)) {
                printf(" initialize_wrap : wrap_name = %s ALREADY EXISTS.",
                       wrap_name);
                ABORT
            }
            wrp = wrp->next;
        }
/*-----------------------------------------------------------------------
!  add a new node to the end of the list
!-----------------------------------------------------------------------*/
        wrp = wrap_data->wrap_head;
        while (wrp != NULL)
            wrp = wrp->next;
        wrp->next = (wrap_node*) tmalloc(sizeof(wrap_node));
        wrp = wrp->next;
        strcpy(wrp->wrap_name, wrap_name);
        wrp->l_wrap = 1;
        wrp->proc_total = 0;
        wrp->send_message_total = 0;
        wrp->recv_message_total = 0;
        wrp->proc = NULL;
        wrp->next = NULL;
    }

    // comm = parallel_get_communicator(communicator_name)
    // make comm = MPI_COMM_WORLD

/*-----------------------------------------------------------------------
!  copy information from the component node to the wrap node
!-----------------------------------------------------------------------*/
    MPI_Comm_size(MPI_COMM_WORLD, &wrp->npe_comm);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrp->rnk_comm);
    wrp->comm      = MPI_COMM_WORLD;
    wrp->flt_type  = MPI_DOUBLE;

    sbdmn_len = grid_params->nsdm;
    wrp->sbdmn_len = sbdmn_len;
    wrp->sbdmn_lst = subdomain_lst;  /* should we make a copy of subdomain_lst? */
    wrp->sbdmn_iota = subdomain_iota;
/*-----------------------------------------------------------------------
!  allocate memory for the sbdmn_perimeter_type and set sbdmn_perimeter_type
!-----------------------------------------------------------------------*/
    
    wrp->sbdmn_perimeter_type   = calloc_2D_int(sbdmn_len, 4);

    wrp->l_sbdmn_pentagon_north = calloc_1D_int(sbdmn_len);
    wrp->l_sbdmn_pentagon_south = calloc_1D_int(sbdmn_len);
    wrp->l_sbdmn_north_pole     = calloc_1D_int(sbdmn_len);
    wrp->l_sbdmn_south_pole     = calloc_1D_int(sbdmn_len);

    set_sbdmn_perimeter_type(wrp, grid_params->grid, grid_params->level_max,
                             level, subdomain_iota, subdomain_lst, sbdmn_len);
/*-----------------------------------------------------------------------
!  set receive lists
!-----------------------------------------------------------------------*/
    set_recv_lst(wrp, grid_params->grid, grid_params->level_max, level,
                 subdomain_iota, subdomain_lst, sbdmn_len);

    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->recv_total > 0) {
            if (proc->proc_nmbr != wrp->rnk_comm)
               wrp->recv_message_total++;

            proc->recv_lst = (int*) tmalloc(4 * proc->recv_total * sizeof(int));
            proc->i1   = proc->recv_lst + proc->recv_total;
            proc->j1   = proc->i1       + proc->recv_total;
            proc->nsd1 = proc->j1       + proc->recv_total;

            instr = proc->recv;
            m = 0;
            while (instr != NULL) {
               proc->recv_lst[m] = instr->tag_glbl;
               proc->i1      [m] = instr->i1;
               proc->j1      [m] = instr->j1;
               proc->nsd1    [m] = instr->nsd1;
               m++;
               instr = instr->next;
            }
        }
        proc = proc->next;
    }
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  COMMUNICATE RECEIVE LISTS (FANCY VERSION)  (091112)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
    wrap_data->msgtag++;

    // the total number of processes the local process communicates with
    proc_total = wrp->proc_total;

    if (proc_total > 0) {
        int *send_total_lst, *recv_total_lst;
        MPI_Request *send_req, *recv_req;
        MPI_Status  *send_status, *recv_status;

        send_total_lst = (int*) malloc(2*proc_total * sizeof(int));
        recv_total_lst = send_total_lst + proc_total;

        send_req = (MPI_Request*) malloc(2*proc_total * sizeof(MPI_Request));
        recv_req = send_req + proc_total;
        for (m=0; m<2*proc_total; m++)
            send_req[m] = MPI_REQUEST_NULL;

        send_status = (MPI_Status*) calloc(2*proc_total, sizeof(MPI_Status));
        recv_status = send_status + proc_total;

#define USE_ISEND_IRECV
#ifdef USE_ISEND_IRECV
/*-----------------------------------------------------------------------
!  set the buffers to receive the number of pieces of information that
!  the local process will send to the neighboring process
!-----------------------------------------------------------------------*/
        proc = wrp->proc;
        for (n=0; n<proc_total; n++) {
            send_total_lst[n] = 0;
            MPI_Irecv(send_total_lst+n, 1, MPI_INT,
                      proc->proc_nmbr, 111, wrp->comm, recv_req+n);
            proc = proc->next;
        }
/*-----------------------------------------------------------------------
!  send the number of pieces of information that the local process
!  expects to receive from the neighboring process
!-----------------------------------------------------------------------*/
        proc = wrp->proc;
        for (n=0; n<proc_total; n++) {
            recv_total_lst[n] = proc->recv_total;
            MPI_Isend(recv_total_lst+n, 1, MPI_INT,
                      proc->proc_nmbr, 111, wrp->comm, send_req+n);
            proc = proc->next;
        }
/*-----------------------------------------------------------------------
!  wait
!-----------------------------------------------------------------------*/
        MPI_Waitall(2*proc_total, send_req, send_status);
#else
        proc = wrp->proc;
        for (n=0; n<proc_total; n++) {
            recv_total_lst[proc->proc_nmbr] = proc->recv_total;
            proc = proc->next;
        }
        MPI_Alltoall(recv_total_lst, 1, MPI_INT, send_total_lst, 1, MPI_INT, wrp->comm);
#endif
/*-----------------------------------------------------------------------
!  allocate memory for the buffer for the global tags to be sent
!-----------------------------------------------------------------------*/
        proc = wrp->proc;
        for (n=0; n<proc_total; n++) {
            proc->recv_msgtag = wrap_data->msgtag;
            proc->send_msgtag = wrap_data->msgtag;
#ifdef USE_ISEND_IRECV
            proc->send_total  = send_total_lst[n];
            proc->send_lst = (int*) tmalloc(send_total_lst[n] * sizeof(int));
#else
            proc->send_total  = send_total_lst[proc->proc_nmbr];
            proc->send_lst = (int*) tmalloc(send_total_lst[proc->proc_nmbr] * sizeof(int));
#endif
            proc = proc->next;
        }
/*-----------------------------------------------------------------------
!  post receives for the global tags to be sent
!-----------------------------------------------------------------------*/
        proc = wrp->proc;
        for (n=0; n<proc_total; n++) {
            MPI_Irecv(proc->send_lst, proc->send_total, MPI_INT,
                      proc->proc_nmbr, 113, wrp->comm, recv_req+n);
            proc = proc->next;
        }
/*-----------------------------------------------------------------------
!  send the list of global tags to be received by the local process
!-----------------------------------------------------------------------*/
        proc = wrp->proc;
        for (n=0; n<proc_total; n++) {
            MPI_Isend(proc->recv_lst, proc->recv_total, MPI_INT,
                      proc->proc_nmbr, 113, wrp->comm, send_req+n);
            proc = proc->next;
        }
/*-----------------------------------------------------------------------
!  wait
!-----------------------------------------------------------------------*/
        MPI_Waitall(2*proc_total, send_req, send_status);
        free(send_status);
        free(send_req);
        free(send_total_lst);
    }    // (proc_total > 0)
    MPI_Barrier(wrp->comm);

/*-----------------------------------------------------------------------
!  set send lists
!-----------------------------------------------------------------------*/
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->send_total > 0) {
            if (proc->proc_nmbr != wrp->rnk_comm)
                wrp->send_message_total++;

            proc->i0 = (int*) tmalloc(3*proc->send_total * sizeof(int));
            proc->j0   = proc->i0 + proc->send_total;
            proc->nsd0 = proc->j0 + proc->send_total;
            for (n=0; n<proc->send_total; n++) {
                get_index(grid_params->grid, grid_params->level_max, level,
                          subdomain_iota, subdomain_lst, sbdmn_len,
                          proc->send_lst[n], ix);  /* ix[] is 1-based */
                proc->i0[n]   = ix[0];
                proc->j0[n]   = ix[1];
                proc->nsd0[n] = ix[2];
            }
        }
        proc = proc->next;
    }
/*-----------------------------------------------------------------------
!  report
!-----------------------------------------------------------------------*/
    if (l_report) {
        FILE *fptr;
        char filename[256];
        sprintf(filename, "./wrp_report/%s_%d_%d_%d", wrap_name, level,
                subdomain_iota, wrp->rnk_comm);
        fptr = fopen(filename, "w");
        fprintf(fptr, "%6d%6d%3d%3d%3d\n",wrp->npe_comm, wrp->rnk_comm,
                wrp->proc_total, wrp->recv_message_total,
                wrp->send_message_total);

        proc_total = wrp->proc_total;
        if (proc_total > 0) {
            proc = wrp->proc;
            for (n=0; n<proc_total; n++) {
                fprintf(fptr,"\n");
                fprintf(fptr,"%6d\n", proc->proc_nmbr);
                fprintf(fptr,"\n");
                fprintf(fptr,"%4d%6d\n", proc->recv_total,proc->recv_msgtag);
                if (proc->recv_total > 0) {
                    for (m=0; m<proc->recv_total; m++)
                        fprintf(fptr,"%6d%4d%4d%3d\n", proc->recv_lst[m],
                                proc->i1[m], proc->j1[m], proc->nsd1[m]);
                }
                fprintf(fptr,"\n");
                fprintf(fptr,"%4d%6d\n", proc->send_total,proc->send_msgtag);
                if (proc->send_total > 0) {
                    for (m=0; m<proc->send_total; m++)
                        fprintf(fptr,"%6d%4d%4d%3d\n", proc->send_lst[m],
                                proc->i0[m], proc->j0[m], proc->nsd0[m]);
                }
                proc = proc->next;
            }
        }
        fclose(fptr);
    }
}

/*----< set_sbdmn_perimeter_type() >------------------------------------------*/
static
void set_sbdmn_perimeter_type(wrap_node *wrp,
                              grid_node *grid,
                              int        level_max,
                              int        level,
                              int        iota,
                              int       *subdomain_lst, /* [sbdmn_len] */
                              int        sbdmn_len)
{
    int i, j, nsd, n, tag_glbl_1, tag_glbl_2;
    grid_node *ptr1, *ptr2, *ptr_nghbr;

    for (i=0; i<sbdmn_len; i++)
        for (j=0; j<4; j++)
            wrp->sbdmn_perimeter_type[i][j] = 0;

    if (level > 0) {

        for (nsd=0; nsd<sbdmn_len; nsd++) {

            /* the values in subdomain_lst[] are 0-based, tags are 1-based */
            tag_glbl_1 = 3+(POWER2(2*(level_max-iota)))*(subdomain_lst[nsd]);
            ptr1 = set_ptr(grid, level_max, level, tag_glbl_1); // bottom-left cell

// subdomain is adjacent to the north pole
            ptr2 = ptr1;
            for (n=1; n<POWER2(level-iota); n++)
                ptr2 = ptr2->nghbr[3];
            wrp->l_sbdmn_north_pole[nsd] = ptr2->nghbr[3]->l_pole_north;

// subdomain is adjacent to the south pole
            ptr2 = ptr1;
            for (n=1; n<POWER2(level-iota); n++)
                ptr2 = ptr2->nghbr[1];
            wrp->l_sbdmn_south_pole[nsd] = ptr2->nghbr[1]->l_pole_south;

            ptr2 = ptr1;
            for (n=1; n<POWER2(level-iota); n++)
                ptr2 = ptr2->nghbr[2];

            tag_glbl_2 = ptr2->tag_glbl; // top-right cell

            wrp->l_sbdmn_pentagon_north[nsd] = ptr1->l_pentagon_north;
            wrp->l_sbdmn_pentagon_south[nsd] = ptr1->l_pentagon_south;

// edge 1 : bottom edge
            if (ptr1->l_pentagon_south)
                wrp->sbdmn_perimeter_type[nsd][0] = 1;
            else {
                ptr_nghbr = ptr1->nghbr[6];
                if (ptr_nghbr->nghbr[3]->tag_glbl != ptr1->tag_glbl)
                    wrp->sbdmn_perimeter_type[nsd][0] = 1;
            }

// edge 2 : right edge
            ptr_nghbr = ptr2->nghbr[1];
            if (ptr_nghbr->nghbr[4]->tag_glbl != ptr2->tag_glbl)
                wrp->sbdmn_perimeter_type[nsd][1] = 2;

// edge 3 : top edge
            ptr_nghbr = ptr2->nghbr[3];
            if (ptr_nghbr->nghbr[6]->tag_glbl != ptr2->tag_glbl)
               wrp->sbdmn_perimeter_type[nsd][2] = 3;

// edge 4 : left edge
            if (ptr1->l_pentagon_north)
                wrp->sbdmn_perimeter_type[nsd][3] = 4;
            else {
                ptr_nghbr = ptr1->nghbr[4];
                if (ptr_nghbr->nghbr[1]->tag_glbl != ptr1->tag_glbl)
                    wrp->sbdmn_perimeter_type[nsd][3] = 4;
            }
        }
    }
}

/*----< set_recv_lst() >------------------------------------------------------*/
static
void set_recv_lst(wrap_node *wrp,
                  grid_node *grid,
                  int        level_max,
                  int        level,
                  int        iota,
                  int       *subdomain_lst,
                  int        sbdmn_len)
{
    int tsklen, el, nsd, tag_glbl, m, n, offset;
    int sbdmn_north[5], sbdmn_south[5], tag_north[5],tag_south[5];
    int *r1, *r2, *n1, *i1, *j1, *di1, *dj1;
    grid_node *ptr0, *ptr1;
    wrap_proc_node *proc;
    wrap_list_node *instr;

    el = POWER2(level-iota);

    tsklen = 6;

    r1  = (int*) malloc(7 * tsklen * sizeof(int));
    r2  = r1  + tsklen;
    n1  = r2  + tsklen;
    i1  = n1  + tsklen;
    j1  = i1  + tsklen;
    di1 = j1  + tsklen;
    dj1 = di1 + tsklen;

    n1[0] = 1; n1[1] = el; n1[2] = el; n1[3] =  1; n1[4] = el; n1[5] = el;
    r1[0] = 0; r1[1] =  1; r1[2] =  3; r1[3] =  0; r1[4] =  4; r1[5] =  6;
    r2[0] = 5; r2[1] =  6; r2[2] =  1; r2[3] =  2; r2[4] =  3; r2[5] =  4;

     i1[0] = 1;  i1[1] = 2;  i1[2] = el+2;  i1[3] = el+2;  i1[4] = el+1;  i1[5] =     1;
     j1[0] = 1;  j1[1] = 1;  j1[2] =    2;  j1[3] = el+2;  j1[4] = el+2;  j1[5] =  el+1;
    di1[0] = 0; di1[1] = 1; di1[2] =    0; di1[3] =    0; di1[4] =   -1; di1[5] =     0;
    dj1[0] = 0; dj1[1] = 0; dj1[2] =    1; dj1[3] =    0; dj1[4] =    0; dj1[5] =    -1;

    for (nsd=0; nsd<sbdmn_len; nsd++) {

        /* the values in subdomain_lst[] are 0-based */
        tag_glbl = 3+(POWER2(2*(level_max-iota)))*(subdomain_lst[nsd]);
        ptr0 = set_ptr(grid, level_max, level, tag_glbl);

        for (m=0; m<tsklen; m++) {
            for (n=0; n<n1[m]; n++) {
                if (ptr0->nghbr[r2[m]] != NULL) {
                    ptr1  = ptr0->nghbr[r2[m]];
                    proc  = get_wrap_proc(wrp, ptr1->proc);
                    instr = get_wrap_list("recv", proc);
                    instr->tag_glbl = ptr1->tag_glbl;
                    instr->i1 = i1[m]+di1[m]*n;
                    instr->j1 = j1[m]+dj1[m]*n;
                    instr->nsd1 = nsd+1;
                }
                if (n+1 < n1[m]) ptr0 = ptr0->nghbr[r1[m]];
            }
        }

/*-----------------------------------------------------------------------
!  north pole
!-----------------------------------------------------------------------*/
        for (m=0; m<5; m++)  /* 0-based */
            sbdmn_north[m] = POWER2(2*iota)*(2*m+1) - 1;

        for (m=0; m<5; m++) /* the values in subdomain_lst[] are 0-based */
            if (sbdmn_north[m] == subdomain_lst[nsd])
                break;
        if (m < 5) {
            ptr0 = set_ptr(grid, level_max, level, 1);
            for (n=0; n<5; n++)
                tag_north[n] = ptr0->nghbr[n+1]->tag_glbl;

            for (n=0; n<5; n++) {
                if (subdomain_lst[nsd] == sbdmn_north[0])
                    break;
                cshift_1D(sbdmn_north, 5, sizeof(int), 1);
                cshift_1D(tag_north,   5, sizeof(int), 1);
            }

            ptr1  = set_ptr(grid, level_max, level, tag_north[2]);
            proc  = get_wrap_proc(wrp, ptr1->proc);
            instr = get_wrap_list("recv", proc);
            instr->tag_glbl = ptr1->tag_glbl;
            instr->i1 = el+2;
            instr->j1 = 1;
            instr->nsd1 = nsd+1;

            ptr1  = set_ptr(grid, level_max, level, tag_north[3]);
            proc  = get_wrap_proc(wrp, ptr1->proc);
            instr = get_wrap_list("recv", proc);
            instr->tag_glbl = ptr1->tag_glbl;
            instr->i1 = 1;
            instr->j1 = el+2;
            instr->nsd1 = nsd+1;
        }
/*-----------------------------------------------------------------------
!  south pole
!-----------------------------------------------------------------------*/
        offset = ((POWER2(iota)+1)*(POWER2(iota)-1))/3;
        for (m=0; m<5; m++)
            sbdmn_south[m] = offset + POWER2(2*iota)*(2*m+1);
        for (m=0; m<5; m++) /* the values in subdomain_lst[] are 0-based */
            if (sbdmn_south[m] == subdomain_lst[nsd])
                break;
        if (m < 5) {
            ptr0 = set_ptr(grid, level_max, level, 2);
            for (n=0; n<5; n++)
                tag_south[n] = ptr0->nghbr[5-n]->tag_glbl;

            for (n=0; n<5; n++) {
                if (subdomain_lst[nsd] == sbdmn_south[0])
                    break;
                cshift_1D(sbdmn_south, 5, sizeof(int), 1);
                cshift_1D(tag_south,   5, sizeof(int), 1);
            }

            ptr1  = set_ptr(grid, level_max, level, tag_south[3]);
            proc  = get_wrap_proc(wrp, ptr1->proc);
            instr = get_wrap_list("recv", proc);
            instr->tag_glbl = ptr1->tag_glbl;
            instr->i1 = el+2;
            instr->j1 = 1;
            instr->nsd1 = nsd+1;

            ptr1  = set_ptr(grid, level_max, level, tag_south[2]);
            proc  = get_wrap_proc(wrp, ptr1->proc);
            instr = get_wrap_list("recv", proc);
            instr->tag_glbl = ptr1->tag_glbl;
            instr->i1 = 1;
            instr->j1 = el+2;
            instr->nsd1 = nsd+1;
        }
    }
    free(r1);
}

/*----< get_wrap_proc() >-----------------------------------------------------*/
static
wrap_proc_node* get_wrap_proc(wrap_node *wrp,
                              int        proc_nmbr)
{
    int l_found;
    wrap_proc_node *proc;

    if (proc_nmbr == RNK_NONEXISTENT) {
        printf(" get_wrap_proc :: proc_nmbr==rnk_nonexistent\n");
        ABORT
    }

    if (wrp->proc == NULL) {
        wrp->proc = (wrap_proc_node*) tmalloc(sizeof(wrap_proc_node));
        wrp->proc_total = 1;
        proc = wrp->proc;
        proc->proc_nmbr   = proc_nmbr;
        proc->send_total  = 0;
        proc->recv_total  = 0;
        proc->send_msgtag = -1;
        proc->recv_msgtag = -1;
        proc->send        = NULL;
        proc->recv        = NULL;
        proc->next        = NULL;
        proc->send_lst    = NULL;
        proc->recv_lst    = NULL;
        proc->i0          = NULL;
    }
    else {
        l_found = 0;
        proc = wrp->proc;
        while (proc != NULL) {
            if (proc->proc_nmbr == proc_nmbr) {
                l_found = 1;
                break;
            }
            proc = proc->next;
        }
        if (! l_found) {
            proc = wrp->proc;
            while (proc->next != NULL)
                proc = proc->next;
            proc->next = (wrap_proc_node*) tmalloc(sizeof(wrap_proc_node));
            wrp->proc_total++;
            proc = proc->next;
            proc->proc_nmbr   = proc_nmbr;
            proc->send_total  = 0;
            proc->recv_total  = 0;
            proc->send_msgtag = -1;
            proc->recv_msgtag = -1;
            proc->send        = NULL;
            proc->recv        = NULL;
            proc->next        = NULL;
            proc->send_lst    = NULL;
            proc->recv_lst    = NULL;
            proc->i0          = NULL;
        }
    }
    return proc;
}

/*----< get_wrap_list() >-----------------------------------------------------*/
static
wrap_list_node* get_wrap_list(char           *option,
                              wrap_proc_node *proc)
{
    wrap_list_node *instr;

    if (!strcmp(option, "send")) {
        proc->send_total++;
        if (proc->send == NULL) {
            proc->send = (wrap_list_node*) tmalloc(sizeof(wrap_list_node));
            instr = proc->send;
            instr->next = NULL;
        }
        else {
            instr = proc->send;
            while (instr->next != NULL)
               instr = instr->next;
            instr->next = (wrap_list_node*) tmalloc(sizeof(wrap_list_node));
            instr = instr->next;
            instr->next = NULL;
        }
    }

    if (!strcmp(option, "recv")) {
        proc->recv_total++;
        if (proc->recv == NULL) {
            proc->recv = (wrap_list_node*) tmalloc(sizeof(wrap_list_node));
            instr = proc->recv;
            instr->next = NULL;
        }
        else {
            instr = proc->recv;
            while (instr->next != NULL)
               instr = instr->next;
            instr->next = (wrap_list_node*) tmalloc(sizeof(wrap_list_node));
            instr = instr->next;
            instr->next = NULL;
        }
    }
    return instr;
}

/*----< find_wrap_node() >----------------------------------------------------*/
static
wrap_node* find_wrap_node(MODULE_wrap_data *wrap_data,
                          const char       *wrap_name)
{
    int l_found;
    wrap_node *wrp=NULL;
/*-----------------------------------------------------------------------
!  find the wrap node associated with wrap_name
!-----------------------------------------------------------------------*/
    if (wrap_data->l_allocate_wrap_head)
        printf("wrap : no wrap nodes\n");
    else {
        wrp = wrap_data->wrap_head;
        l_found = 0;
        while (wrp != NULL) {
            if (!strcmp(wrp->wrap_name, wrap_name)) {
                l_found = 1;
                break;
            }
            wrp = wrp->next;
        }
        if (! l_found) {
            printf("wrap: the wrap instructions = %s has not been initialized.",
                   wrap_name);
            ABORT
        }
    }
/*-----------------------------------------------------------------------
!  if SCM then return
!-----------------------------------------------------------------------*/
    if (! wrp->l_wrap) return NULL;

    return wrp;
}

/*----< wrap_face() >---------------------------------------------------------*/
void wrap_face(MODULE_wrap_data    *wrap_data,
               const char          *wrap_name,
               double           ****face,
               const int           *dimlen)  /* [4] */
{
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    adjust_face(wrp, face, dimlen);
    wrap_4d(wrp, face, dimlen);
}

/*----< wrap_vrtx() >---------------------------------------------------------*/
void wrap_vrtx(MODULE_wrap_data      *wrap_data,
               const char            *wrap_name,
               double           ******vrtx,
               const int             *dimlen)  /* [6] */
{
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    wrap_6d(wrp, vrtx, dimlen);
    adjust_vrtx(wrp, vrtx, dimlen);
}

/*----< wrap_edge() >---------------------------------------------------------*/
void wrap_edge(MODULE_wrap_data      *wrap_data,
               const char            *wrap_name,
               double           ******edge,
               const int             *dimlen)  /* [6] */
{
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    wrap_6d(wrp, edge, dimlen);
    adjust_edge(wrp, edge, dimlen);
}

/*----< wrap_face_1lyr() >----------------------------------------------------*/
void wrap_face_1lyr(MODULE_wrap_data    *wrap_data,
                    const char          *wrap_name,
                    double           ****face_1lyr, /* process only [][][][idx] */
                    const int           *dimlen3,   /* [3] */
                    int                  idx)       /* the index of 4th dim */
{
    int i, j, k, dimlen4[4];
    double ****temp_rk4, ***face_1lyr3D;
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    dimlen4[0] = dimlen3[0];
    dimlen4[1] = 1;
    dimlen4[2] = dimlen3[1];
    dimlen4[3] = dimlen3[2];
    temp_rk4 = calloc_4D_dbl(dimlen4[0], dimlen4[1], dimlen4[2], dimlen4[3]);

    if (idx < 0) {
        face_1lyr3D = *face_1lyr;
        for (k=0; k<dimlen3[0]; k++)
            for (j=0; j<dimlen3[1]; j++)
                for (i=0; i<dimlen3[2]; i++)
                    temp_rk4[k][0][j][i] = face_1lyr3D[k][j][i];
    }
    else {
        for (k=0; k<dimlen3[0]; k++)
            for (j=0; j<dimlen3[1]; j++)
                for (i=0; i<dimlen3[2]; i++)
                    temp_rk4[k][0][j][i] = face_1lyr[k][j][i][idx];
    }

    adjust_face(wrp, temp_rk4, dimlen4);
    wrap_4d(wrp, temp_rk4, dimlen4);

    if (idx < 0) {
        for (k=0; k<dimlen3[0]; k++)
            for (j=0; j<dimlen3[1]; j++)
                for (i=0; i<dimlen3[2]; i++)
                    face_1lyr3D[k][j][i] = temp_rk4[k][0][j][i];
    } else {
        for (k=0; k<dimlen3[0]; k++)
            for (j=0; j<dimlen3[1]; j++)
                for (i=0; i<dimlen3[2]; i++)
                    face_1lyr[k][j][i][idx] = temp_rk4[k][0][j][i];
    }

    free_4D_dbl(temp_rk4);
}

/*----< wrap_vrtx_1lyr() >----------------------------------------------------*/
void wrap_vrtx_1lyr(MODULE_wrap_data     *wrap_data,
                    const char           *wrap_name,
                    double           *****vrtx_1lyr,
                    const int            *dimlen5) /* [5] */
{
    int g, h, i, j, k, dimlen6[6];
    double ******temp_rk6;
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    dimlen6[0] = dimlen5[0];
    dimlen6[1] = 1;
    dimlen6[2] = dimlen5[1];
    dimlen6[3] = dimlen5[2];
    dimlen6[4] = dimlen5[3];
    dimlen6[5] = dimlen5[4];
    temp_rk6 = calloc_6D_dbl(dimlen6[0], dimlen6[1], dimlen6[2],
                             dimlen6[3], dimlen6[4], dimlen6[5]);
    for (k=0; k<dimlen5[0]; k++)
        for (j=0; j<dimlen5[1]; j++)
            for (i=0; i<dimlen5[2]; i++)
                for (h=0; h<dimlen5[3]; h++)
                    for (g=0; g<dimlen5[4]; g++)
                        temp_rk6[k][0][j][i][h][g] = vrtx_1lyr[k][j][i][h][g];

    wrap_6d(wrp, temp_rk6, dimlen6);
    adjust_vrtx(wrp, temp_rk6, dimlen6);

    for (k=0; k<dimlen5[0]; k++)
        for (j=0; j<dimlen5[1]; j++)
            for (i=0; i<dimlen5[2]; i++)
                for (h=0; h<dimlen5[3]; h++)
                    for (g=0; g<dimlen5[4]; g++)
                        vrtx_1lyr[k][j][i][h][g] = temp_rk6[k][0][j][i][h][g];

    free_6D_dbl(temp_rk6);
}

/*----< wrap_vrtx_scalar() >--------------------------------------------------*/
void wrap_vrtx_scalar(MODULE_wrap_data     *wrap_data,
                      const char           *wrap_name,
                      double           *****vrtx_scalar,
                      const int            *dimlen5) /* [5] */
{
    int g, h, i, j, k, dimlen6[6];
    double ******temp_rk6;
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    dimlen6[0] = dimlen5[0];
    dimlen6[1] = dimlen5[1];
    dimlen6[2] = dimlen5[2];
    dimlen6[3] = dimlen5[3];
    dimlen6[4] = dimlen5[4];
    dimlen6[5] = 1;
    temp_rk6 = calloc_6D_dbl(dimlen6[0], dimlen6[1], dimlen6[2],
                             dimlen6[3], dimlen6[4], dimlen6[5]);
    for (k=0; k<dimlen5[0]; k++)
        for (j=0; j<dimlen5[1]; j++)
            for (i=0; i<dimlen5[2]; i++)
                for (h=0; h<dimlen5[3]; h++)
                    for (g=0; g<dimlen5[4]; g++)
                        temp_rk6[k][j][i][h][g][0] = vrtx_scalar[k][j][i][h][g];

    wrap_6d(wrp, temp_rk6, dimlen6);
    adjust_vrtx(wrp, temp_rk6, dimlen6);

    for (k=0; k<dimlen5[0]; k++)
        for (j=0; j<dimlen5[1]; j++)
            for (i=0; i<dimlen5[2]; i++)
                for (h=0; h<dimlen5[3]; h++)
                    for (g=0; g<dimlen5[4]; g++)
                        vrtx_scalar[k][j][i][h][g] = temp_rk6[k][j][i][h][g][0];

    free_6D_dbl(temp_rk6);
}

/*----< wrap_vrtx_scalar_1lyr() >---------------------------------------------*/
void wrap_vrtx_scalar_1lyr(MODULE_wrap_data    *wrap_data,
                           const char          *wrap_name,
                           double           ****vrtx_scalar_1lyr,
                           const int           *dimlen4) /* [4] */
{
    int h, i, j, k, dimlen6[6];
    double ******temp_rk6;
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    dimlen6[0] = dimlen4[0];
    dimlen6[1] = 1;
    dimlen6[2] = dimlen4[1];
    dimlen6[3] = dimlen4[2];
    dimlen6[4] = dimlen4[3];
    dimlen6[5] = 1;
    temp_rk6 = calloc_6D_dbl(dimlen6[0], dimlen6[1], dimlen6[2],
                             dimlen6[3], dimlen6[4], dimlen6[5]);
    for (k=0; k<dimlen4[0]; k++)
        for (j=0; j<dimlen4[1]; j++)
            for (i=0; i<dimlen4[2]; i++)
                for (h=0; h<dimlen4[3]; h++)
                    temp_rk6[k][0][j][i][h][0] = vrtx_scalar_1lyr[k][j][i][h];

    wrap_6d(wrp, temp_rk6, dimlen6);
    adjust_vrtx(wrp, temp_rk6, dimlen6);

    for (k=0; k<dimlen4[0]; k++)
        for (j=0; j<dimlen4[1]; j++)
            for (i=0; i<dimlen4[2]; i++)
                for (h=0; h<dimlen4[3]; h++)
                    vrtx_scalar_1lyr[k][j][i][h] = temp_rk6[k][0][j][i][h][0];

    free_6D_dbl(temp_rk6);
}

/*----< wrap_edge_1lyr() >----------------------------------------------------*/
void wrap_edge_1lyr(MODULE_wrap_data      *wrap_data,
                    const char            *wrap_name,
                    double            *****edge_1lyr,
                    const int             *dimlen5) /* [5] */
{
    int g, h, i, j, k, dimlen6[6];
    double ******temp_rk6;
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    dimlen6[0] = dimlen5[0];
    dimlen6[1] = 1;
    dimlen6[2] = dimlen5[1];
    dimlen6[3] = dimlen5[2];
    dimlen6[4] = dimlen5[3];
    dimlen6[5] = dimlen5[4];
    temp_rk6 = calloc_6D_dbl(dimlen6[0], dimlen6[1], dimlen6[2],
                             dimlen6[3], dimlen6[4], dimlen6[5]);
    for (k=0; k<dimlen5[0]; k++)
        for (j=0; j<dimlen5[1]; j++)
            for (i=0; i<dimlen5[2]; i++)
                for (h=0; h<dimlen5[3]; h++)
                    for (g=0; g<dimlen5[4]; g++)
                        temp_rk6[k][0][j][i][h][g] = edge_1lyr[k][j][i][h][g];

    wrap_6d(wrp, temp_rk6, dimlen6);
    adjust_edge(wrp, temp_rk6, dimlen6);

    for (k=0; k<dimlen5[0]; k++)
        for (j=0; j<dimlen5[1]; j++)
            for (i=0; i<dimlen5[2]; i++)
                for (h=0; h<dimlen5[3]; h++)
                    for (g=0; g<dimlen5[4]; g++)
                        edge_1lyr[k][j][i][h][g] = temp_rk6[k][0][j][i][h][g];

    free_6D_dbl(temp_rk6);
}

/*----< wrap_edge_scalar() >--------------------------------------------------*/
void wrap_edge_scalar(MODULE_wrap_data     *wrap_data,
                      const char           *wrap_name,
                      double           *****edge_scalar,
                      const int            *dimlen5)  /* [5] */
{
    int g, h, i, j, k, dimlen6[6];
    double ******temp_rk6;
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    dimlen6[0] = dimlen5[0];
    dimlen6[1] = dimlen5[1];
    dimlen6[2] = dimlen5[2];
    dimlen6[3] = dimlen5[3];
    dimlen6[4] = dimlen5[4];
    dimlen6[5] = 1;
    temp_rk6 = calloc_6D_dbl(dimlen6[0], dimlen6[1], dimlen6[2],
                             dimlen6[3], dimlen6[4], dimlen6[5]);
    for (k=0; k<dimlen5[0]; k++)
        for (j=0; j<dimlen5[1]; j++)
            for (i=0; i<dimlen5[2]; i++)
                for (h=0; h<dimlen5[3]; h++)
                    for (g=0; g<dimlen5[4]; g++)
                        temp_rk6[k][j][i][h][g][0] = edge_scalar[k][j][i][h][g];

    wrap_6d(wrp, temp_rk6, dimlen6);
    adjust_edge(wrp, temp_rk6, dimlen6);

    for (k=0; k<dimlen5[0]; k++)
        for (j=0; j<dimlen5[1]; j++)
            for (i=0; i<dimlen5[2]; i++)
                for (h=0; h<dimlen5[3]; h++)
                    for (g=0; g<dimlen5[4]; g++)
                        edge_scalar[k][j][i][h][g] = temp_rk6[k][j][i][h][g][0];

    free_6D_dbl(temp_rk6);
}

/*----< wrap_edge_scalar_1lyr() >---------------------------------------------*/
void wrap_edge_scalar_1lyr(MODULE_wrap_data    *wrap_data,
                           const char          *wrap_name,
                           double           ****edge_scalar_1lyr,
                           const int           *dimlen4) /* [4] */
{
    int h, i, j, k, dimlen6[6];
    double ******temp_rk6;
    wrap_node *wrp;

    wrp = find_wrap_node(wrap_data, wrap_name);
    if (wrp == NULL) return;

    dimlen6[0] = dimlen4[0];
    dimlen6[1] = 1;
    dimlen6[2] = dimlen4[1];
    dimlen6[3] = dimlen4[2];
    dimlen6[4] = dimlen4[3];
    dimlen6[5] = 1;
    temp_rk6 = calloc_6D_dbl(dimlen6[0], dimlen6[1], dimlen6[2],
                             dimlen6[3], dimlen6[4], dimlen6[5]);
    for (k=0; k<dimlen4[0]; k++)
        for (j=0; j<dimlen4[1]; j++)
            for (i=0; i<dimlen4[2]; i++)
                for (h=0; h<dimlen4[3]; h++)
                    temp_rk6[k][0][j][i][h][0] = edge_scalar_1lyr[k][j][i][h];

    wrap_6d(wrp, temp_rk6, dimlen6);
    adjust_edge(wrp, temp_rk6, dimlen6);

    for (k=0; k<dimlen4[0]; k++)
        for (j=0; j<dimlen4[1]; j++)
            for (i=0; i<dimlen4[2]; i++)
                for (h=0; h<dimlen4[3]; h++)
                    edge_scalar_1lyr[k][j][i][h] = temp_rk6[k][0][j][i][h][0];

    free_6D_dbl(temp_rk6);
}


/*----< wrap_4d() >-----------------------------------------------------------*/
static
void wrap_4d(wrap_node    *wrp,
             double    ****x1,
             const int    *dimlen)  /* [4] */
{
    int i, nlyr, nsdm, count, m, n, iota, k, nsd;
    int sbdmn_pentagon_north[5], sbdmn_pentagon_south[5];
    MPI_Request *req;
    MPI_Status *status;
    wrap_proc_node *proc;

    nsdm = dimlen[0];
    nlyr = dimlen[1];

    int total_messages = wrp->send_message_total + wrp->recv_message_total;
    req = (MPI_Request*) malloc(total_messages * sizeof(MPI_Request));
    for (i=0; i<total_messages; i++)
        req[i] = MPI_REQUEST_NULL;
    status = (MPI_Status*) calloc(total_messages, sizeof(MPI_Status));

/*-----------------------------------------------------------------------
!  receives
!-----------------------------------------------------------------------*/
    // if (l_wrptmr) timer (event_name="wrp_recv",action="start")
    count = 0;
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->recv_total > 0 && proc->proc_nmbr != wrp->rnk_comm) {
            proc->recv_rk2 = calloc_2D_dbl(nlyr, proc->recv_total);
            MPI_Irecv(proc->recv_rk2[0], proc->recv_total*nlyr,
                      wrp->flt_type, proc->proc_nmbr, proc->recv_msgtag,
                      wrp->comm, &req[count++]);
        }
        proc = proc->next;
    }
    // IF (l_wrptmr) CALL timer (event_name="wrp_recv",action="stop")
/*-----------------------------------------------------------------------
!  sends
!-----------------------------------------------------------------------*/
    // IF (l_wrptmr) CALL timer (event_name="wrp_send",action="start")
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->send_total > 0 && proc->proc_nmbr != wrp->rnk_comm) {
            proc->send_rk2 = calloc_2D_dbl(nlyr, proc->send_total);
            for (k=0; k<nlyr; k++)
                for (m=0; m<proc->send_total; m++)
                   proc->send_rk2[k][m] =
                   x1[proc->nsd0[m]-1][k][proc->j0[m]-1][proc->i0[m]-1];
            MPI_Isend(proc->send_rk2[0], proc->send_total*nlyr,
                      wrp->flt_type, proc->proc_nmbr, proc->send_msgtag,
                      wrp->comm, &req[count++]);
        }
        proc = proc->next;
    }
    // IF (l_wrptmr) CALL timer (event_name="wrp_send",action="stop")
/*-----------------------------------------------------------------------
!  wrap local
!-----------------------------------------------------------------------*/
    // IF (l_wrptmr) CALL timer (event_name="wrp_locl",action="start")
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->proc_nmbr == wrp->rnk_comm) {
            for (k=0; k<nlyr; k++)
                for (m=0; m<proc->recv_total; m++)
                    x1[proc->nsd1[m]-1][k][proc->j1[m]-1][proc->i1[m]-1] =
                    x1[proc->nsd0[m]-1][k][proc->j0[m]-1][proc->i0[m]-1];
         }
         proc = proc->next;
    }
    // IF (l_wrptmr) CALL timer (event_name="wrp_locl",action="stop")
/*-----------------------------------------------------------------------
!  wait
!-----------------------------------------------------------------------*/
    // IF (l_wrptmr) CALL timer (event_name="wrp_wait",action="start")

    MPI_Waitall(count, req, status);

    // IF (l_wrptmr) CALL timer (event_name="wrp_wait",action="stop")
/*-----------------------------------------------------------------------
!  unpack
!-----------------------------------------------------------------------*/
    // IF (l_wrptmr) CALL timer (event_name="wrp_unpk",action="start")
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->recv_total > 0) {
            if (proc->proc_nmbr != wrp->rnk_comm) {
                for (k=0; k<nlyr; k++)
                    for (m=0; m<proc->recv_total; m++)
                        x1[proc->nsd1[m]-1][k][proc->j1[m]-1][proc->i1[m]-1] =
                        proc->recv_rk2[k][m];
            }
        }
        proc = proc->next;
    }
    // IF (l_wrptmr) CALL timer (event_name="wrp_unpk",action="stop")
/*-----------------------------------------------------------------------
!  deallocate memory
!-----------------------------------------------------------------------*/
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->send_total > 0 && proc->proc_nmbr != wrp->rnk_comm)
            free_2D_dbl(proc->send_rk2);
        if (proc->recv_total > 0 && proc->proc_nmbr != wrp->rnk_comm)
            free_2D_dbl(proc->recv_rk2);
        proc = proc->next;
    }
/*-----------------------------------------------------------------------
!  special
!-----------------------------------------------------------------------*/
    iota = wrp->sbdmn_iota;
    for (n=0; n<5; n++)
        sbdmn_pentagon_north[n] = POWER2(2*iota)*(2*n);

    for (n=0; n<5; n++)
        sbdmn_pentagon_south[n] = POWER2(2*iota)*(2*n+1);

    /* wrp->sbdmn_lst is an array of size gcrm_param->nsdm, values are 0-based */
    for (nsd=0; nsd<wrp->sbdmn_len; nsd++) {
        for (n=0; n<5; n++) {
            if (sbdmn_pentagon_north[n] == wrp->sbdmn_lst[nsd]) {
                for (k=0; k<nlyr; k++)
                    x1[nsd][k][1][0] = x1[nsd][k][0][0];
                break;
            }
         }
        for (n=0; n<5; n++) {
            if (sbdmn_pentagon_south[n] == wrp->sbdmn_lst[nsd]) {
                for (k=0; k<nlyr; k++)
                    x1[nsd][k][0][1] = x1[nsd][k][0][0];
                break;
            }
         }
    }

    free(req);
    free(status);
}

/*----< wrap_6d() >-----------------------------------------------------------*/
static
void wrap_6d(wrap_node      *wrp,
             double    ******x2,
             const int      *dimlen)  /* [6] */
{
    int n1, n2, nlyr, nsdm, count, m, n, iota, k, nsd;
    int sbdmn_pentagon_north[5], sbdmn_pentagon_south[5];
    MPI_Request *req;
    MPI_Status *status;
    wrap_proc_node *proc;

    n1   = dimlen[5];
    n2   = dimlen[4];
    nlyr = dimlen[1];
    nsdm = dimlen[0];

    int total_messages = wrp->send_message_total + wrp->recv_message_total;
    req = (MPI_Request*) malloc(total_messages * sizeof(MPI_Request));
    for (n=0; n<total_messages; n++)
        req[n] = MPI_REQUEST_NULL;
    status = (MPI_Status*) calloc(total_messages, sizeof(MPI_Status));

/*-----------------------------------------------------------------------
!  receives
!-----------------------------------------------------------------------*/
    count = 0;
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->recv_total > 0 && proc->proc_nmbr != wrp->rnk_comm) {
            proc->recv_rk4 = calloc_4D_dbl(nlyr, proc->recv_total, n2, n1);
            MPI_Irecv(proc->recv_rk4[0][0][0], n1*n2*proc->recv_total*nlyr,
                      wrp->flt_type, proc->proc_nmbr, proc->recv_msgtag,
                      wrp->comm, &req[count++]);
        }
        proc = proc->next;
    }
/*-----------------------------------------------------------------------
!  sends
!-----------------------------------------------------------------------*/
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->send_total > 0 && proc->proc_nmbr != wrp->rnk_comm) {
            proc->send_rk4 = calloc_4D_dbl(nlyr, proc->send_total, n2, n1);
            for (k=0; k<nlyr; k++)
                for (m=0; m<proc->send_total; m++)
                    memcpy(proc->send_rk4[k][m][0], 
                           x2[proc->nsd0[m]-1][k][proc->j0[m]-1][proc->i0[m]-1][0],
                           n1*n2*sizeof(double));
            MPI_Isend(proc->send_rk4[0][0][0], n1*n2*proc->send_total*nlyr,
                      wrp->flt_type, proc->proc_nmbr, proc->send_msgtag,
                      wrp->comm, &req[count++]);
        }
        proc = proc->next;
    }
/*-----------------------------------------------------------------------
!  wrap local
!-----------------------------------------------------------------------*/
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->proc_nmbr == wrp->rnk_comm) {
            for (k=0; k<nlyr; k++)
                for (m=0; m<proc->recv_total; m++)
                    if (proc->nsd1[m] != proc->nsd0[m] || proc->j1[m] != proc->j0[m] || proc->i1[m] != proc->i0[m])
                        memcpy(x2[proc->nsd1[m]-1][k][proc->j1[m]-1][proc->i1[m]-1][0],
                               x2[proc->nsd0[m]-1][k][proc->j0[m]-1][proc->i0[m]-1][0],
                               n1*n2*sizeof(double));
         }
         proc = proc->next;
    }
/*-----------------------------------------------------------------------
!  wait
!-----------------------------------------------------------------------*/
    MPI_Waitall(total_messages, req, status);
/*-----------------------------------------------------------------------
!  unpack
!-----------------------------------------------------------------------*/
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->recv_total > 0) {
            if (proc->proc_nmbr != wrp->rnk_comm) {
                for (k=0; k<nlyr; k++)
                    for (m=0; m<proc->recv_total; m++)
                        memcpy(x2[proc->nsd1[m]-1][k][proc->j1[m]-1][proc->i1[m]-1][0],
                               proc->recv_rk4[k][m][0], n1*n2*sizeof(double));
            }
        }
        proc = proc->next;
    }
/*-----------------------------------------------------------------------
!  deallocate memory
!-----------------------------------------------------------------------*/
    proc = wrp->proc;
    while (proc != NULL) {
        if (proc->send_total > 0 && proc->proc_nmbr != wrp->rnk_comm)
            free_4D_dbl(proc->send_rk4);
        if (proc->recv_total > 0 && proc->proc_nmbr != wrp->rnk_comm)
            free_4D_dbl(proc->recv_rk4);
        proc = proc->next;
    }
/*-----------------------------------------------------------------------
!  special
!-----------------------------------------------------------------------*/
    iota = wrp->sbdmn_iota;
    for (n=0; n<5; n++)
        sbdmn_pentagon_north[n] = POWER2(2*iota)*(2*n);

    for (n=0; n<5; n++)
        sbdmn_pentagon_south[n] = POWER2(2*iota)*(2*n+1);

    for (nsd=0; nsd<wrp->sbdmn_len; nsd++) {
        for (n=0; n<5; n++) { /* values of sbdmn_lst[] are 0-based */
            if (sbdmn_pentagon_north[n] == wrp->sbdmn_lst[nsd]) {
                for (k=0; k<nlyr; k++)
                    memcpy(x2[nsd][k][1][0][0],
                           x2[nsd][k][0][0][0], n1*n2*sizeof(double));
                break;
            }
         }
        for (n=0; n<5; n++) { /* values of sbdmn_lst[] are 0-based */
            if (sbdmn_pentagon_south[n] == wrp->sbdmn_lst[nsd]) {
                for (k=0; k<nlyr; k++)
                    memcpy(x2[nsd][k][0][1][0],
                           x2[nsd][k][0][0][0], n1*n2*sizeof(double));
                break;
            }
         }
    }

    free(req);
    free(status);
}

/*----< adjust_face() >-------------------------------------------------------*/
static
void adjust_face(wrap_node    *wrp,
                 double    ****face,
                 const int    *dimlen) /* [4] */
{
    int z, y;
    for (z=0; z<dimlen[0]; z++)
        for (y=0; y<dimlen[1]; y++) {
            face[z][y][dimlen[2]-1][0] = face[z][y][1][1];
            face[z][y][0][dimlen[3]-1] = face[z][y][1][1];
        }
}


/*----< adjust_edge() >-------------------------------------------------------*/
static
void adjust_edge(wrap_node       *wrp,
                 double    ******edge,
                 const int      *dimlen) /* [6] */
{
    int im, jm, km, lm, i, j, k, l, nsd;
/*
!   INTEGER (KIND=int_kind) :: &
!      iota,n,sbdmn_pentagon_north(5),sbdmn_pentagon_south(5)
!.:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:. .:.

!   iota = wrp->sbdmn_iota
!   sbdmn_pentagon_north(:)=(/ (1+2**(2*iota)*n,n=0,8,2) /)
!   sbdmn_pentagon_south(:)=(/ (1+2**(2*iota)*n,n=1,9,2) /)

   im = SIZE (edge,DIM=3); jm = SIZE (edge,DIM=4);
*/
    im = dimlen[3];
    jm = dimlen[2];
    km = dimlen[1];
    lm = dimlen[5];

    for (nsd=0; nsd<wrp->sbdmn_len; nsd++) {
        /* BOTTOM EDGE */
        if (wrp->sbdmn_perimeter_type[nsd][0] == 1) {
            if (wrp->l_sbdmn_pentagon_south[nsd]) {
                for (i=1; i<im-1; i++) {
                    for (k=0; k<km; k++)
                        for (l=0; l<lm; l++)
                            edge[nsd][k][0][i-1][0][l] = edge[nsd][k][0][i][2][l];

                    for (k=0; k<km; k++)
                         cshift_down_2D(&edge[nsd][k][0][i][0][0],
                                        3, lm*sizeof(double));
                }
                for (k=0; k<km; k++)
                    for (l=0; l<lm; l++)
                        edge[nsd][k][0][0][0][l] = 0.0;
            }
            else {
                for (k=0; k<km; k++)
                    cshift_down_2D(&edge[nsd][k][0][0][0][0],
                                   3, lm*sizeof(double));
                for (i=1; i<im-1; i++) {
                    for (k=0; k<km; k++)
                        for (l=0; l<lm; l++)
                            edge[nsd][k][0][i-1][0][l] = edge[nsd][k][0][i][2][l];

                    for (k=0; k<km; k++)
                         cshift_down_2D(&edge[nsd][k][0][i][0][0],
                                        3, lm*sizeof(double));
                }
            }
        }
        /* RIGHT EDGE */
        if (wrp->sbdmn_perimeter_type[nsd][1] == 2) {
            for (k=0; k<km; k++)
                cshift_up_2D(&edge[nsd][k][1][im-1][0][0],
                             3, lm*sizeof(double));
            for (j=2; j<jm; j++) {
                for (k=0; k<km; k++)
                    for (l=0; l<lm; l++)
                        edge[nsd][k][j][im-1][2][l] = edge[nsd][k][j][im-1][0][l];

                for (k=0; k<km; k++)
                    cshift_up_2D(&edge[nsd][k][j][im-1][0][0],
                                 3, lm*sizeof(double));
            }
        }
        /* TOP EDGE */
        if (wrp->sbdmn_perimeter_type[nsd][2] == 3) {
            for (k=0; k<km; k++)
                cshift_down_2D(&edge[nsd][k][jm-1][1][0][0],
                               3, lm*sizeof(double));

            for (i=2; i<im; i++) {
                for (k=0; k<km; k++)
                    for (l=0; l<lm; l++)
                        edge[nsd][k][jm-1][i-2][0][l] = edge[nsd][k][jm-1][i][2][l];

                for (k=0; k<km; k++)
                    cshift_down_2D(&edge[nsd][k][jm-1][i][0][0],
                                   3, lm*sizeof(double));
            }
        }
        /* LEFT EDGE */
        if (wrp->sbdmn_perimeter_type[nsd][3] == 4) {
            if (wrp->l_sbdmn_pentagon_north[nsd]) {
                for (j=2; j<jm; j++) {
                    for (k=0; k<km; k++)
                        for (l=0; l<lm; l++)
                            edge[nsd][k][j-1][0][2][l] = edge[nsd][k][j][0][0][l];
    
                    for (k=0; k<km; k++)
                        cshift_up_2D(&edge[nsd][k][j][0][0][0],
                                     3, lm*sizeof(double));
                }
                for (k=0; k<km; k++)
                    for (l=0; l<lm; l++)
                        edge[nsd][k][0][0][2][l] = 0.0;
            }
            else {
                for (k=0; k<km; k++)
                    cshift_up_2D(&edge[nsd][k][0][0][0][0],
                                 3, lm*sizeof(double));
                for (j=2; j<jm; j++) {
                    for (k=0; k<km; k++)
                        for (l=0; l<lm; l++)
                            edge[nsd][k][j-1][0][2][l] = edge[nsd][k][j][0][0][l];
    
                    for (k=0; k<km; k++)
                        cshift_up_2D(&edge[nsd][k][j][0][0][0],
                                     3, lm*sizeof(double));
                }
            }
        }
        /* NORTH AND SOUTH POLE */
        /* north pole */
        for (k=0; k<km; k++)
            cshift_down_2D(&edge[nsd][k][jm-1][0][0][0],
                           3, lm*sizeof(double));
        /* south pole */
        for (k=0; k<km; k++)
            cshift_up_2D(&edge[nsd][k][0][im-1][0][0],
                         3, lm*sizeof(double));
    }
}

/*----< adjust_vrtx() >-------------------------------------------------------*/
static
void adjust_vrtx(wrap_node      *wrp,
                 double    ******vrtx,
                 const int      *dimlen) /* [6] */
{
    int im, jm, i, j, k, nsd;

    im = dimlen[3];
    jm = dimlen[2];

    for (nsd=0; nsd<wrp->sbdmn_len; nsd++) {
/*-----------------------------------------------------------------------
!  BOTTOM EDGE
!-----------------------------------------------------------------------*/
        if (wrp->sbdmn_perimeter_type[nsd][0] == 1) {
            if (wrp->l_sbdmn_pentagon_south[nsd]) {
                for (j=0; j<dimlen[1]; j++) {
                    for (i=0; i<dimlen[5]; i++) {
                        // vrtx(:,2,1,1,:,nsd) = vrtx(:,2,1,1,:,nsd)
                        vrtx[nsd][j][0][1][1][i] = vrtx[nsd][j][0][0][0][i];
                    }
                }
                for (i=2; i<im-1; i++) {
                    for (k=0; k<dimlen[1]; k++) {
                        for (j=0; j<dimlen[5]; j++) {
                            vrtx[nsd][k][0][i-1][0][j] = vrtx[nsd][k][0][i][1][j];
                            vrtx[nsd][k][0][i  ][1][j] = vrtx[nsd][k][0][i][0][j];
                        }
                    }
                }
            }
            else {
                for (i=0; i<im-2; i++) {
                    for (k=0; k<dimlen[1]; k++) {
                        for (j=0; j<dimlen[5]; j++) {
                            vrtx[nsd][k][0][i][1][j] = vrtx[nsd][k][0][i  ][0][j];
                            vrtx[nsd][k][0][i][0][j] = vrtx[nsd][k][0][i+1][1][j];
                        }
                    }
                }
                for (k=0; k<dimlen[1]; k++) {
                    for (j=0; j<dimlen[5]; j++) {
                        vrtx[nsd][k][0][im-2][1][j] = vrtx[nsd][k][0][im-2][0][j];
                    }
                }
            }
        }
/*-----------------------------------------------------------------------
!  RIGHT EDGE
!-----------------------------------------------------------------------*/
        if (wrp->sbdmn_perimeter_type[nsd][1] == 2) {
            for (j=1; j<jm-2; j++) {
                for (k=0; k<dimlen[1]; k++) {
                    for (i=0; i<dimlen[5]; i++) {
                        vrtx[nsd][k][j][im-1][0][i] = vrtx[nsd][k][j  ][im-1][1][i];
                        vrtx[nsd][k][j][im-1][1][i] = vrtx[nsd][k][j+1][im-1][0][i];
                    }
                }
            }
        }
/*-----------------------------------------------------------------------
!  TOP EDGE
!-----------------------------------------------------------------------*/
        if (wrp->sbdmn_perimeter_type[nsd][2] == 3) {
            for (i=1; i<im-2; i++) {
                for (k=0; k<dimlen[1]; k++) {
                    for (j=0; j<dimlen[5]; j++) {
                        vrtx[nsd][k][jm-1][i][1][j] = vrtx[nsd][k][jm-1][i  ][0][j];
                        vrtx[nsd][k][jm-1][i][0][j] = vrtx[nsd][k][jm-1][i+1][1][j];
                    }
                }
            }
        }
/*-----------------------------------------------------------------------
!  LEFT EDGE
!-----------------------------------------------------------------------*/
        if (wrp->sbdmn_perimeter_type[nsd][3] == 4) {
            if (wrp->l_sbdmn_pentagon_north[nsd]) {
                // vrtx(:,1,1,1,:,nsd) = vrtx(:,1,1,1,:,nsd)
                for (k=0; k<dimlen[1]; k++)
                    for (j=0; j<dimlen[5]; j++)
                        vrtx[nsd][k][1][0][0][j] = vrtx[nsd][k][0][0][1][j];

                for (j=1; j<jm-2; j++) {
                    for (k=0; k<dimlen[1]; k++) {
                        for (i=0; i<dimlen[5]; i++) {
                            vrtx[nsd][k][j-1][0][1][i] = vrtx[nsd][k][j][0][0][i];
                            vrtx[nsd][k][j  ][0][0][i] = vrtx[nsd][k][j][0][1][i];
                        }
                    }
                }
            }
            else {
                for (j=1; j<jm-2; j++) {
                    for (k=0; k<dimlen[1]; k++) {
                        for (i=0; i<dimlen[5]; i++) {
                            vrtx[nsd][k][j][0][0][i] = vrtx[nsd][k][j  ][0][1][i];
                            vrtx[nsd][k][j][0][1][i] = vrtx[nsd][k][j+1][0][0][i];
                        }
                    }
                }
                for (k=0; k<dimlen[1]; k++) {
                    for (i=0; i<dimlen[5]; i++) {
                        vrtx[nsd][k][jm-2][0][0][i] = vrtx[nsd][k][jm-2][0][1][i];
                    }
                }
            }
        }
/*-----------------------------------------------------------------------
!  NORTH POLE. see wrap_convention_corner_north_pole.png
!-----------------------------------------------------------------------*/
        if (wrp->l_sbdmn_north_pole[nsd]) {
            for (k=0; k<dimlen[1]; k++) {
                for (i=0; i<dimlen[5]; i++) {
                    vrtx[nsd][k][jm-2][0][1][i] = vrtx[nsd][k][jm-1][0   ][1][i];
                    vrtx[nsd][k][jm-1][0][0][i] = vrtx[nsd][k][0   ][im-1][1][i];
                    vrtx[nsd][k][jm-1][0][1][i] = vrtx[nsd][k][0   ][im-1][0][i];
                    vrtx[nsd][k][jm-1][1][1][i] = -999.0;  /* undefined */
                }
            }
        }
/*-----------------------------------------------------------------------
!  SOUTH POLE. see wrap_convention_corner_south_pole.png
!-----------------------------------------------------------------------*/
        if (wrp->l_sbdmn_south_pole[nsd]) {
            for (k=0; k<dimlen[1]; k++) {
                for (i=0; i<dimlen[5]; i++) {
                    vrtx[nsd][k][0][im-2][0][i] = vrtx[nsd][k][0   ][im-1][0][i];
                    vrtx[nsd][k][0][im-1][0][i] = vrtx[nsd][k][jm-1][0   ][1][i];
                    vrtx[nsd][k][0][im-1][1][i] = vrtx[nsd][k][jm-1][0   ][0][i];
                    vrtx[nsd][k][1][im-1][0][i] = -999.0;  /* undefined */
                }
            }
        }
    }
}

