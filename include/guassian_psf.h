//
// Created by peo on 17-6-27.
//

#ifndef PSF_GUASSIAN_PSF_H
#define PSF_GUASSIAN_PSF_H
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdio>

#ifndef M_PI
#define M_PI (3.1415926535897932384626433832795)
#endif
#define M_2PI (6.283185307179586476925286766559)

#define BESSEL_LEN 1010 /* length of Bessel function lookup table */
#define BESSEL_RES 1.0 /* resolution of Bessel function lookup table */
#define BESSEL_INT 60		/* steps for integrating the Bessel integral */

/* lookup table for Bessel function of the first kind, orders 0, 1, and 2 */
double bessel_lut[BESSEL_LEN*3];

/*****************************************************************************/
/* C functions
 */

#define floor_int(x) (int)floor((double)(x))
#define ceil_int(x) (int)ceil((double)(x))


/* Return Bessel function x for orders 0, 1, and 2.
 * Values are linearly interpolated from the lookup table.
 */
void bessel_lookup(double x, double* result);

/* Initialize the global Bessel function lookup table.
 * Values are approximated by intergrating Bessel's integral.
 */
int bessel_init(void);

/* Caculate 2D Gaussian distribution.
 */
int gaussian2d(double *out, int *shape, double* sigma);

/* Caculate Gaussian parameters for the nonparaxial widefield case.
 */
void sigma_widefield(double* sz, double* sr, double nk, double cosa);

/* Calculate Gaussian parameters for points spread function approximation
 * according to B Zhang et al. Appl. Optics (46) 1819-29, 2007.
 */
int gaussian_sigma(double* sz, double* sr, double lex, double lem,
                   double NA, double n, double r, int widefield, int paraxial);
/* Apodization function for exitation.
 */
double apodization_excitation(double ct, double st, double _, double beta);

/* Apodization function for isotropic fluorescence emission.
 */
double apodization_emission(double ct, double st, double M, double _);

/* Caculate the Point Spread Function for unpolarized or circular polarized
 * light according
 */
int psf(
        int type,        /* PSF type: 0: excitation, 1: emission */
        double *data,    /* output array[shape[0]][shape[1]] */
        int* shape,      /* shape of data array */
        double* uvdim,   /* optical units in u and v dimension */
        double M,        /* lateral magnification factor */
        double sinalpha, /* numerical aperture / refractive index of medium */
        double beta,     /* underfilling ratio (1.0) */
        double gamma,    /* ex_wavelen / em_wavelen / refractive index (1.0) */
        int intsteps     /* number of steps for integrating over theta (50) */
);


/* Calculate the observation volume for unpolarized light by multiplying the
 * excitation PSF with the convolution of the emission PSF and detector kernel.
 *
 * The PSFs and the detector kernel must have equal physical sizes per pixel.
 * The PSFs must be in zr space, the detector kernel in xy space.
 */
int obsvol(
        int dimz,         /* size of PSF arrays in z dimension */
        int dimr,         /* size of PSF arrays in r dimension */
        int dimd,         /* pixel size of detector array (must be square) */
        double *obsvol,   /* output array[dimu, dimr, dimr] */
        double *ex_psf,   /* excitation PSF array[dimu, dimr] */
        double *em_psf,   /* emission PSF array[dimu, dimr] */
        double *detector); /* detector kernel array[dimd, dimd] or NULL */

/*
 *Calculate the detector kernel for integration over pinhole with trapezoid rule.
 *
 *The radius denotes the outer radius, except for square shape, where it denotes
 *the inner radius.
 */
int pinhole_kernel(int corners, double* out, int dim, double radius);
/* Apply rotational symmetry around the 1st dimension.
 */
int zr2zxy(double* data, double* out, int dimz, int dimr);
void psf_generate();

#endif //PSF_GUASSIAN_PSF_H
