/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_utilities.c 1858 2013-06-12 04:49:45Z wkliao $
 */

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>  /* memset() */
#include <math.h>

#include "grid_utilities.h"
#include "util.h"

#define PI 3.14159265358979323846264338327950288419716939937510

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif


/*----< mid_point() >---------------------------------------------------------*/
void mid_point(double *p1,     /* [3] */
               double *p2,     /* [3] */
               double *p_out)  /* [3] */
{
/*
! PURPOSE : given two points (p1,p2) measured in (x,y,z) and lying on
!           the unit sphere, find the mid-point of the arch between p1 
!           and p2.
!
! NOTE : p1 and p2 are assumed to be of unit length
! NOTE : p_out is normalized to unit length
*/
    p_out[0] = 0.5 * (p1[0] + p2[0]);
    p_out[1] = 0.5 * (p1[1] + p2[1]);
    p_out[2] = 0.5 * (p1[2] + p2[2]);
    unit_vector(p_out, p_out);
}

/*----< tangent_to_sphere() >-------------------------------------------------*/
void tangent_to_sphere(double *p1,        /* [3] */
                       double *p2,        /* [3] */
                       double *p_out)     /* [3] */
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : given two points (p1,p2) measured in (x,y,z) and lying on
!           the unit sphere, find the vector (p_out) that lies 
!           in the plane tangent to the unit sphere at the point p1 
!           and points in the direction of the projection of p2 onto 
!           the tangent plane.
!
! NOTE : p1 and p2 are assumed to be of unit length
! NOTE : p_out is normalized to unit length
!.......................................................................*/

    double tmp1[3], tmp2[3];
    cross_product(p1, p2, tmp1);
    cross_product(tmp1, p1, tmp2);
    unit_vector(tmp2, p_out);
}

/*----< voronoi_corner() >----------------------------------------------------*/
void voronoi_corner(double *p1,       /* [3] */
                    double *p2,       /* [3] */
                    double *p3,       /* [3] */
                    double *p_out)    /* [3] */
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : given three points (p1,p2,p3) measured in (x,y,z) and lying
!           on the unit sphere, find the point (p_out) that is equidistant
!           from p1, p2 and p3.
!
! NOTE : p1, p2 and p3 are assumed to be given in a counterclockwise fashion
!.......................................................................*/
    double temp1[3], temp2[3], temp3[3];

    temp1[0] = p2[0] - p1[0];
    temp1[1] = p2[1] - p1[1];
    temp1[2] = p2[2] - p1[2];
    temp2[0] = p3[0] - p1[0];
    temp2[1] = p3[1] - p1[1];
    temp2[2] = p3[2] - p1[2];

    cross_product(temp1, temp2, temp3);
    unit_vector(temp3, p_out);
}

/*----< cross_product() >-----------------------------------------------------*/
void cross_product(double *p1,        /* [3] */
                   double *p2,        /* [3] */
                   double *p_out)     /* [3] */
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : given (p1,p2) measured in (x,y,z) and lying on the unit sphere, 
!           compute their cross product.  result returned in p_out
!.......................................................................*/

    if ((p_out <= p1 - 3 || p1 + 3 <= p_out) &&
        (p_out <= p2 - 3 || p2 + 3 <= p_out)) {
        /* if p_out is neither p1 nor p2 */
        p_out[0] = p1[1]*p2[2] - p1[2]*p2[1];
        p_out[1] = p1[2]*p2[0] - p1[0]*p2[2];
        p_out[2] = p1[0]*p2[1] - p1[1]*p2[0];
    }
    else {
        double temp[3];
        temp[0] = p1[1]*p2[2] - p1[2]*p2[1];
        temp[1] = p1[2]*p2[0] - p1[0]*p2[2];
        temp[2] = p1[0]*p2[1] - p1[1]*p2[0];

        p_out[0] = temp[0];
        p_out[1] = temp[1];
        p_out[2] = temp[2];
    }
}

/*----< unit_vector() >-------------------------------------------------------*/
void unit_vector(double *p,        /* [3] */
                 double *p_out)    /* [3] */
{
    double dist = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];

    if (dist == 0.0) {
        /* this must be the case p[] is a zero vector */
        if (p != p_out)
            p_out[0] = p_out[1] = p_out[2] = 0.0;
        return;
    }
    dist = sqrt(dist);
    p_out[0] = p[0] / dist;
    p_out[1] = p[1] / dist;
    p_out[2] = p[2] / dist;
}

/*----< lonlat_to_xyz() >-----------------------------------------------------*/
void lonlat_to_xyz(double *lonlat,   /* [2] */
                   double *p_out)    /* [3] */
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : given a point measured in longitude and latitude, 
!           find the Cartesian coordinates p = (x,y,z).
!
! NOTE : x,y,z are in the range [-1,1].
!
! NOTE : the two coordinate systems are related in the following manner:
!          (0   , 0    )  (greenwich -- equator) transforms to (1, 0, 0)
!          (PI/2, 0    )  (90¼ east  -- equator) transforms to (0, 1, 0) 
!          (0   , PI/2 )  (north pole)           transforms to (0, 0, 1)
!.......................................................................*/

    p_out[0] = cos(lonlat[0]) * cos(lonlat[1]);
    p_out[1] = sin(lonlat[0]) * cos(lonlat[1]);
    p_out[2] = sin(lonlat[1]);
}


/*----< spherical_triangle_area() >-------------------------------------------*/
double spherical_triangle_area(double *p1,  /* in: [3] */
                               double *p2,  /* in: [3] */
                               double *p3)  /* in: [3] */
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : given three points measured in (x,y,z) lying on the unit
!           sphere, find the spherical area of the triangle formed by
!           connecting these three points
!
! NOTE : p1, p2, p3 are assumed to be of unit length
! NOTE : p_out has units of radians^2
!
! NOTE : see Mathematical CRC Tables for formula
!        area = angle1 + angle2 + angle3 - PI
!        where angle1, angle2 are the vertex angles
!.......................................................................*/
    double t1[3], t2[3], dot, angle1, angle2, angle3;

    tangent_to_sphere(p1, p2, t1);
    tangent_to_sphere(p1, p3, t2);
    dot    = MAX(-1.0, MIN(1.0, t1[0]*t2[0]+t1[1]*t2[1]+t1[2]*t2[2]));
    angle1 = acos(dot);

    tangent_to_sphere(p2, p3, t1);
    tangent_to_sphere(p2, p1, t2);
    dot    = MAX(-1.0, MIN(1.0, t1[0]*t2[0]+t1[1]*t2[1]+t1[2]*t2[2]));
    angle2 = acos(dot);

    tangent_to_sphere(p3, p1, t1);
    tangent_to_sphere(p3, p2, t2);
    dot    = MAX(-1.0, MIN(1.0, t1[0]*t2[0]+t1[1]*t2[1]+t1[2]*t2[2]));
    angle3 = acos(dot);

    return (angle1 + angle2 + angle3 - PI);
}

/*----< area_corner_kites() >-------------------------------------------------*/
void area_corner_kites(double *p1,      /* [3] */
                       double *p2,      /* [3] */
                       double *p3,      /* [3] */
                       double *area_c)  /* [3] out */
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : given three points p1, p2 and p3 measured in (x,y,z) lying 
!           on the unit sphere, find the spherical area formed by
!           connecting these three points.  the total area is returned
!           in three kites. 
!           see /Users/ross/swm/figures/archive/corner_area.eps
!
! NOTE : p1, p2, p3 are assumed to be of unit length
! NOTE : p1, p2, p3 are assumed to be given in a counterclockwise fashion
! NOTE : p_out has units of radians^2
!
! NOTE : see Mathematical CRC Tables for formula
!        area = angle1 + angle2 + angle3 - PI
!        where angle1, angle2 are the vertex angles
!.......................................................................*/
    double mid12[3], mid23[3], mid31[3], mid123[3];

    mid_point (p1, p2, mid12);
    mid_point (p2, p3, mid23);
    mid_point (p3, p1, mid31);

    voronoi_corner(p1, p2, p3, mid123);

    area_c[0] = spherical_triangle_area(    p1, mid12, mid31) +
                spherical_triangle_area(mid123, mid31, mid12);
    area_c[1] = spherical_triangle_area(    p2, mid23, mid12) +
                spherical_triangle_area(mid123, mid12, mid23);
    area_c[2] = spherical_triangle_area(    p3, mid31, mid23) +
                spherical_triangle_area(mid123, mid23, mid31);
}

/*----< arch_distance() >-----------------------------------------------------*/
double arch_distance(double *p1,     /* [3] */
                     double *p2)     /* [3] */
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : given two points measured in (x,y,z) lying on the unit
!           sphere, find the distance measured along the surface of the
!           sphere between the two points.
!!.......................................................................*/
    double dot;

    dot = p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];

    if (dot >= 1.0) dot =  1.0;
    if (dot <=-1.0) dot = -1.0;

    return acos(dot);
}

/*----< set_vctr_wghts_crn() >------------------------------------------------*/
void set_vctr_wghts_crn(double  *p1,       /* [3] */
                        double  *p2,       /* [3] */
                        double  *p3,       /* [3] */
                        double **weights)  /* [2][3] out */
{
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PURPOSE : assume points p1, p2 and p3 (measured in (x,y,z)) are on
!           the unit sphere.  compute the weights associated with the 
!           gradient operator defined at cell corners.  
!           see Ringler and Randall, A Potential Enstrophy and Energy 
!           Conserving Numerical Scheme for Solution of the Shallow-Water 
!           Equations on a Geodesic Grid, Mon. Wea. Rev. 130, 1397-1410, 
!           equation (40) and (41)
!.......................................................................*/
    int i, n;
    double north_pole[3], local_vert[3], crn[3], north[3], east[3];
    double mid[3], tng[3], lngth[3], tmp[3];
    double trans_mat[3][3], p[3][3], nrm[3][3];

    north_pole[0] = 0.0;
    north_pole[1] = 0.0;
    north_pole[2] = 1.0;
    local_vert[0] = 0.0;
    local_vert[1] = 0.0;
    local_vert[2] = 1.0;

    voronoi_corner(p1, p2, p3, crn);
    tangent_to_sphere(crn, north_pole, north);
    cross_product(north, crn, east);
    unit_vector(east, east);

    for (n=0; n<3; n++) {
        trans_mat[n][0] =  east[n];
        trans_mat[n][1] = north[n];
        trans_mat[n][2] =   crn[n];
    }

    memcpy(p[0], p1, 3*sizeof(double));
    memcpy(p[1], p2, 3*sizeof(double));
    memcpy(p[2], p3, 3*sizeof(double));

    for (n=0; n<3; n++) {
        mid_point(p[0], p[1], mid);
        tangent_to_sphere(crn, mid, tng);
        for (i=0; i<3; i++) tng[i] = -tng[i];
        mat_vec_mul(3, 3, trans_mat, tng, tmp);
        unit_vector(tmp, tng);

        cross_product(tng, local_vert, tmp);
        unit_vector(tmp, nrm[n]);
        lngth[n] = arch_distance(crn, mid);

        cshift_up_2D(&p[0][0], 3, 3*sizeof(double));
    }

    for (n=0; n<3; n++) {
        for (i=0; i<2; i++)
            weights[i][n] = lngth[2]*nrm[2][i]+lngth[0]*(-nrm[0][i]);

        cshift_1D(lngth, 3, sizeof(double), 1);
        cshift_up_2D(&nrm[0][0], 3, 3*sizeof(double));
    }
}


/*----< mat_vec_mul() >-------------------------------------------------------*/
/* matrix-vector multiplication */
void mat_vec_mul(int      m,
                 int      n,
                 double   x[m][n],  /* in:  [m][n] */
                 double  *y,        /* in:  [n] */
                 double  *z)        /* out: [m] */
{
    int i, j;
    for (i=0; i<m; i++) {
        z[i] = 0.0;
        for (j=0; j<n; j++) {
            z[i] += x[i][j] * y[j];
        }
    }
}

/*----< matmul() >------------------------------------------------------------*/
/* 2D matrix multiplication */
void matmul(int      m,
            int      q,
            int      n,
            double **x,  /* in:  [m][q] */
            double **y,  /* in:  [q][n] */
            double **z)  /* out: [m][n] */
{
    int i, j, k;
    for (i=0; i<m; i++) {
        for (j=0; j<n; j++) {
            z[i][j] = 0.0;
            for (k=0; k<q; k++)
                z[i][j] += x[i][k] * y[k][j];
        }
    }
}

