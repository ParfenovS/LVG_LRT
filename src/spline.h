/* interpolation/akima.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman
 */
#pragma once
#include <stdlib.h>
#include <math.h>
#include <vector>

class akima_spline // see https://raw.githubusercontent.com/ampl/gsl/master/interpolation/akima.c
{
private:

    size_t size;
    bool allocated;
    double * stateb;
    double * statec;
    double * stated;
    double * state_m;
    double * x_array;
    double * y_array;

    /* evaluation accelerator */
    typedef struct {
        size_t  cache;        /* cache of index   */
        size_t  miss_count;   /* keep statistics  */
        size_t  hit_count;
    } gsl_interp_accel;

    gsl_interp_accel *acc;

    gsl_interp_accel *gsl_interp_accel_alloc()
    {
        gsl_interp_accel *a = (gsl_interp_accel *) malloc (sizeof (gsl_interp_accel));
        a->cache = 0;
        a->hit_count = 0;
        a->miss_count = 0;
        return a;
    }

    size_t gsl_interp_bsearch(double x, size_t index_lo, size_t index_hi)
    {
        size_t ilo = index_lo;
        size_t ihi = index_hi;
        while (ihi > ilo + 1) {
            size_t i = (ihi + ilo)/2;
            if(x_array[i] > x) ihi = i;
            else ilo = i;
        }
        return ilo;
    }

    size_t gsl_interp_accel_find(size_t len, double x)
    {
        size_t x_index = acc->cache;
        if (x < x_array[x_index]) {
            acc->miss_count++;
            acc->cache = gsl_interp_bsearch(x, 0, x_index);
        } else if (x >= x_array[x_index + 1]) {
            acc->miss_count++;
            acc->cache = gsl_interp_bsearch(x, x_index, len-1);
        }
        else {
            acc->hit_count++;
        }
        return acc->cache;
    }

    void akima_calc(double b[],  double c[],  double d[], size_t size, double m[])
    {
        size_t i;

        for (i = 0; i < (size - 1); i++) {
            const double NE = fabs (m[i + 1] - m[i]) + fabs (m[i - 1] - m[i - 2]);
            if (NE == 0.0) {
                b[i] = m[i];
                c[i] = 0.0;
                d[i] = 0.0;
            } else {
                const double h_i = x_array[i + 1] - x_array[i];
                const double NE_next = fabs (m[i + 2] - m[i + 1]) + fabs (m[i] - m[i - 1]);
                const double alpha_i = fabs (m[i - 1] - m[i - 2]) / NE;
                double alpha_ip1;
                double tL_ip1;
                if (NE_next == 0.0) {
                    tL_ip1 = m[i];
                } else {
                    alpha_ip1 = fabs (m[i] - m[i - 1]) / NE_next;
                    tL_ip1 = (1.0 - alpha_ip1) * m[i] + alpha_ip1 * m[i + 1];
                }
                b[i] = (1.0 - alpha_i) * m[i - 1] + alpha_i * m[i];
                c[i] = (3.0 * m[i] - 2.0 * b[i] - tL_ip1) / h_i;
                d[i] = (b[i] + tL_ip1 - 2.0 * m[i]) / (h_i * h_i);
            }
        }
    }

public:

    akima_spline()
    {
        size = 0;
        allocated = false;
    }

    void akima_init(const vector <double> & inp_x_array, const vector <double> & inp_y_array)
    {
        size = inp_x_array.size();
        stateb = (double *) malloc (size * sizeof (double));
        statec = (double *) malloc (size * sizeof (double));
        stated = (double *) malloc (size * sizeof (double));
        state_m = (double *) malloc ((size + 4) * sizeof (double));
        acc = gsl_interp_accel_alloc ();

        x_array = (double *) malloc (size * sizeof (double));
        y_array = (double *) malloc (size * sizeof (double));
        for (size_t iar = 0; iar < size; iar++) {
            x_array[iar] = inp_x_array[iar];
            y_array[iar] = inp_y_array[iar];
        }

        allocated = true;

        double * m = state_m + 2; /* offset so we can address the -1,-2
                                        components */

        size_t i;
        for (i = 0; i <= size - 2; i++) {
            m[i] = (y_array[i + 1] - y_array[i]) / (x_array[i + 1] - x_array[i]);
        }
        
        /* non-periodic boundary conditions */
        m[-2] = 3.0 * m[0] - 2.0 * m[1];
        m[-1] = 2.0 * m[0] - m[1];
        m[size - 1] = 2.0 * m[size - 2] - m[size - 3];
        m[size] = 3.0 * m[size - 2] - 2.0 * m[size - 3];
        
        akima_calc (stateb, statec, stated, size, m);
    }

    void clear()
    {
        if (allocated) {
            free (stateb);
            free (statec);
            free (stated);
            free (state_m);
            free (acc);
            free (x_array);
            free (y_array);
        }
        allocated = false;
    }

    ~akima_spline()
    {
        clear();
    }

    double akima_eval(double x)
    {
        size_t index;
        
        if (acc != 0) {
            index = gsl_interp_accel_find (size, x);
        } else {
            index = gsl_interp_bsearch (x, 0, size - 1);
        }
        
        /* evaluate */
        {
            const double x_lo = x_array[index];
            const double delx = x - x_lo;
            const double b = stateb[index];
            const double c = statec[index];
            const double d = stated[index];
            return y_array[index] + delx * (b + delx * (c + d * delx));
        }
    }

    double akima_eval_deriv(double x)
    {
        size_t index;
        
        if (acc != 0) {
            index = gsl_interp_accel_find (size, x);
        } else {
            index = gsl_interp_bsearch (x, 0, size - 1);
        }
        
        /* evaluate */
        {
            double x_lo = x_array[index];
            double delx = x - x_lo;
            double b = stateb[index];
            double c = statec[index];
            double d = stated[index];
            return b + delx * (2.0 * c + 3.0 * d * delx);
        }
    }
};
