/*
 *  Copyright (C) 2013, Northwestern University.
 *  See COPYRIGHT notice COPYING in top-level directory.
 *
 *  $Id: grid_utilities.h 1858 2013-06-12 04:49:45Z wkliao $
 */

#ifndef H_GRID_UTILITIES
#define H_GRID_UTILITIES

void mid_point(double *p1, double *p2, double *p_out);
void tangent_to_sphere(double *p1, double *p2, double *p_out);
void voronoi_corner(double *p1, double *p2, double *p3, double *p_out);
void cross_product(double *p1, double *p2, double *p_out);
void unit_vector(double *p, double *p_out);
void lonlat_to_xyz(double *lonlat, double *p_out);
double spherical_triangle_area(double *p1, double *p2, double *p3);
void area_corner_kites(double *p1, double *p2, double *p3, double *area_c);
double arch_distance(double *p1, double *p2);
void set_vctr_wghts_crn(double *p1, double *p2, double *p3, double **weights);
void matmul(int m, int q, int n, double **x, double **y, double **z);
void mat_vec_mul(int m, int n, double x[m][n], double *y, double *z);

#endif
