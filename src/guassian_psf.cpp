#include "guassian_psf.h"

/* Return Bessel function x for orders 0, 1, and 2.
 * Values are linearly interpolated from the lookup table.
 */
void bessel_lookup(double x, double* result)
{
	double *p;
	double alpha = x * BESSEL_RES;
	int index = floor_int(alpha);

	if (index < BESSEL_LEN) {
		p = &bessel_lut[index*3];
		alpha -= (double)index;
		result[0] = p[0] + alpha * (p[3] - p[0]);
		result[1] = p[1] + alpha * (p[4] - p[1]);
		result[2] = p[2] + alpha * (p[5] - p[2]);
	} else {
		result[0] = result[1] = result[2] = 0;
	}
}


/* Initialize the global Bessel function lookup table.
 * Values are approximated by intergrating Bessel's integral.
 */
int bessel_init(void)
{
	double x, t, xst, dt, *ptr;
	
	memset(bessel_lut, 0, BESSEL_LEN*3*sizeof(double));

	dt = M_PI / (BESSEL_INT);
	ptr = bessel_lut;
	for (int i=0; i < BESSEL_LEN; i++)
	{
		x = -(double)i / BESSEL_RES;
		for (int j = 0; j < BESSEL_INT; j++)
		{
			t = j * dt;
			xst = x * sin(t);
			ptr[0] += cos(xst);
			ptr[2] += cos(2.0*t + xst);
		}
		ptr[0] /= BESSEL_INT;
		ptr[2] /= BESSEL_INT;
		ptr += 3;
	}

	dt = M_PI / (BESSEL_INT-1);
	ptr = bessel_lut+1;
	for (int i = 0; i < BESSEL_LEN; i++)
	{
		x = (double)i / BESSEL_RES;
		for (int j = 0; j < BESSEL_INT; j++)
		{
			t = j * dt;
			*ptr += cos(t - x*sin(t));
		}
		*ptr /= BESSEL_INT - 1;
		ptr += 3;
	}
	return 0;
}


/* Caculate 2D Gaussian distribution.
 */
int gaussian2d(double *out, int *shape, double* sigma)
{
	double sz, sr, t;
	if ((out == NULL) || (sigma[0] == 0) || (sigma[1] == 0))
	{
		return -1;
	}

	sz = -0.5 / (sigma[0]*sigma[0]);
	sr = -0.5 / (sigma[1]*sigma[1]);

	for (int z = 0; z < shape[0]; z++)
	{
		t = z*z * sz;
		for (int r = 0; r < shape[1]; r++)
		{
			*out++ = exp(t + r*r*sr);
		}
	}
	return 0;
}

/* Caculate Gaussian parameters for the nonparaxial widefield case.
 */
void sigma_widefield(double* sz, double* sr, double nk, double cosa)
{
	double t = pow(cosa, 1.5);
	*sr = 1.0 / (nk * sqrt((4. - 7.*t + 3.*pow(cosa, 3.5)) / (7.*(1. - t))));
	*sz = (5.*sqrt(7.)*(1. - t)) / (nk * sqrt(6.*(4.*pow(cosa, 5) - 
				25.*pow(cosa, 3.5) + 42.*pow(cosa, 2.5) - 25.*t + 4.)));
}

/* Calculate Gaussian parameters for points spread function approximation
 * according to B Zhang et al. Appl. Optics (46) 1819-29, 2007.
 */
int gaussian_sigma(double* sz, double* sr, double lex, double lem,
										double NA, double n, double r, int widefield, int paraxial)
{
	if ((NA <= 0.0) || (n <= 0.0) || (lem <= 0.0) || ((NA/n) >= 1.0))
		return -1;

	if (widefield)
	{
		if (paraxial)
		{
			*sr = sqrt(2.) * lem / (M_2PI * NA);
			*sz = sqrt(6.) * lem / (M_2PI * NA*NA) * n * 2.;
		}
		else
		{
			sigma_widefield(sz, sr, n*M_2PI/lem, cos(asin(NA/n)));
		}
	}
	else
	{
		if ((r <= 0.0) || (lex <= 0.0))
			return -1;
		if(paraxial)
		{
			double kem = M_2PI / lem;
			double c1 = M_2PI / lex * r * NA;
			double c2 = M_2PI / lem * r * NA;
			double J0, J1, J[3];
			bessel_lookup(c2, J);
			J0 = J[0];
			J1 = J[1];
			*sr = sqrt(2./(c1*c1/(r*r) + 
					(4.*c2*J0*J1 - 8.*J1*J1) / (r*r*(J0*J0 + J1*J1 - 1.))));
			*sz = 2.*sqrt(6./((c1*c1*NA*NA)/(r*r*n*n) - 
					(48.*c2*c2*(J0*J0 + J1*J1) - 192.*J1*J1) / 
					(n*n*kem*kem*r*r*r*r*(J0*J0 + J1*J1 - 1.))));
		}
		else
		{
			double e, sz_em, sr_em, sz_ex, sr_ex;
			double cosa = cos(asin(NA/n));
			sigma_widefield(&sz_ex, &sr_ex, n*M_2PI/lex, cosa);
			sigma_widefield(&sz_em, &sr_em, n*M_2PI/lem, cosa);
			e = sr_em*sr_em;
			e = 2.0 * e * e * (exp(r*r/(2.0*e)) - 1.0);
			*sr = sqrt((sr_ex*sr_ex * e) / (e + r*r * sr_ex*sr_ex));
			*sz = sz_ex*sz_em / sqrt(sz_ex*sz_ex + sz_em*sz_em);
		}
	}
	return 0;
}

/* Apodization function for exitation.
 */
double apodization_excitation(double ct, double st, double _, double beta)
{
	return sqrt(ct) * exp(st*st * beta);
}

/* Apodization function for isotropic fluorescence emission.
 */
double apodization_emission(double ct, double st, double M, double _)
{
	double t = M*st;
	return sqrt(ct / sqrt(1.0 - t*t));
}

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
)
{
	int i, j, k, u_shape, v_shape;
	double u, v, t, st, ct, re, im, re0, im0, re1, im1, re2, im2;
    double const0, const1, u_delta, v_delta, bessel[3];
    double alpha; /* integration over theta upper limit */
    double delta; /* step size for integrating over theta */
    double *cache, *cptr, *dptr;
    double (*apodization)(double, double, double, double);

	if ((intsteps < 4) || (sinalpha <= 0.0) || (sinalpha >= 1.0))
        return -1;
	
	switch (type) {
        case 0: /* excitation */
            apodization = apodization_excitation;
            alpha = asin(sinalpha);
            beta = -beta*beta / (sinalpha*sinalpha);
            gamma = M = 1.0;
            break;
        case 1: /* emission */
            apodization = apodization_emission;
            alpha = asin(sinalpha / M);
            beta = 1.0;
            break;
        default:
            return -1;
    }

	delta = alpha / (double)(intsteps-1);

    cache = cptr = (double *)malloc((intsteps*5)*sizeof(double));
    if (cache == NULL)
        return -1;

    const0 = gamma / sinalpha;
    const1 = gamma / (sinalpha * sinalpha);

    /* cache some values used in inner integration loop */
    for (k = 0; k < intsteps; k++) {
        t = k * delta;
        st = sin(t);
        ct = cos(t);
        t = st * apodization(ct, st, M, beta);
        cptr[0] = st * const0;
        cptr[1] = ct * const1;
        cptr[2] = t * (1.0 + ct);
        cptr[3] = t * st * (2.0); /* 4*I1(u,v) */
        cptr[4] = t * (1.0 - ct);
        cptr += 5;
    }

    u_shape = shape[0];
    v_shape = shape[1];
    u_delta = uvdim[0] / (double)(u_shape-1);
    v_delta = uvdim[1] / (double)(v_shape-1);
    dptr = data;

	for (i = 0; i < u_shape; i++) {
        u = u_delta * (double) i;
        for (j = 0; j < v_shape; j++) {
            v = v_delta * (double) j;
            re0 = im0 = re1 = im1 = re2 = im2 = 0.0;
            cptr = cache;
            /* integrate over theta using trapezoid rule */
            bessel_lookup(v * cptr[0], bessel);
            ct = u * cptr[1]; re = cos(ct); im = sin(ct);
            t = bessel[0]*cptr[2]*0.5; re0 += re*t; im0 += im*t;
            t = bessel[1]*cptr[3]*0.5; re1 += re*t; im1 += im*t;
            t = bessel[2]*cptr[4]*0.5; re2 += re*t; im2 += im*t;
            cptr += 5;
            for (k = 1; k < intsteps-1; k++) {
                bessel_lookup(v * cptr[0], bessel);
                ct = u * cptr[1];
                re = cos(ct); /* complex exponential with re=0 */
                im = sin(ct);
                t = bessel[0]*cptr[2]; re0 += re*t; im0 += im*t;
                t = bessel[1]*cptr[3]; re1 += re*t; im1 += im*t;
                t = bessel[2]*cptr[4]; re2 += re*t; im2 += im*t;
                cptr += 5;
            }
            bessel_lookup(v * cptr[0], bessel);
            ct = u * cptr[1]; re = cos(ct); im = sin(ct);
            t = bessel[0]*cptr[2]*0.5; re0 += re*t; im0 += im*t;
            t = bessel[1]*cptr[3]*0.5; re1 += re*t; im1 += im*t;
            t = bessel[2]*cptr[4]*0.5; re2 += re*t; im2 += im*t;

            *dptr++ = (re0*re0 + im0*im0 +
                       re1*re1 + im1*im1 +
                       re2*re2 + im2*im2);
        }
    }
    t = data[0];
    for (i = 0; i < u_shape*v_shape; i++)
        data[i] /= t;

    free(cache);
    return 0;
}


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
    double *detector) /* detector kernel array[dimd, dimd] or NULL */
{
    int z, x, y, r, xx, i, ii, ri, index, indey, *_ri;
    double sum, rf, x2, t;
    double *exptr, *emptr, *ovptr, *_rf, *_em;
    int _dimd = dimd;

    if (detector == NULL) {
        /* approximation for large pinholes == widefield */
        i = 0;
        for (z = 0; z < dimz; z++) {
            sum = em_psf[i] * M_PI * 0.25;
            ii = i;
            for (r = 1; r < dimr; r++) {
                sum += em_psf[++ii] * (double)r;
            }
            sum *= M_2PI;
            for (r = 0; r < dimr; r++, i++) {
                obsvol[i] = ex_psf[i] * sum;
            }
        }
    } else if (dimd < 2) {
        /* approximation for very small pinholes */
        for (i = 0; i < dimz*dimr; i++) {
            obsvol[i] = ex_psf[i] * em_psf[i];
        }
    } else {
        /* use detector/pinhole kernel */
        if (dimd > dimr) dimd = dimr;
        /* floor integer and remainder float of radius at xy */
        _ri = (int *)malloc((dimr*dimd)*sizeof(int));
        if (_ri == NULL)
            return -1;
        _rf = (double *)malloc((dimr*dimd)*sizeof(double));
        if (_rf == NULL) {
            free(_ri);
            return -1;
        }
        /* em_psf at xy */
        _em = (double *)malloc((dimr*dimd)*sizeof(double));
        if (_em == NULL) {
            free(_ri); free(_rf);
            return -1;
        }

        for (x = 0; x < dimd; x++) {
            x2 = (double)(x*x);
            indey = x;
            index = x*dimd;
            _ri[index] = _ri[indey] = x;
            _rf[index] = _rf[indey] = 0.0;
            for (y = 1; y <= x; y++) {
                index++;
                indey += dimd;
                rf = sqrt(x2 + y*y);
                ri = floor_int(rf);
                _ri[index] = _ri[indey] = (dimr > ri) ? ri : -1;
                _rf[index] = _rf[indey] = (dimr > ri+1) ? rf-(double)ri : 0.0;
            }
        }
        for (x = dimd; x < dimr; x++) {
            index = x*dimd;
            _ri[index] = x;
            _rf[index] = 0.0;
            x2 = (double)(x*x);
            for (y = 1; y < dimd; y++) {
                index++;
                rf = sqrt(x2 + y*y);
                ri = floor_int(rf);
                _ri[index] = (dimr > ri) ? ri : -1;
                _rf[index] = (dimr > ri+1) ? rf-(double)ri : 0.0;
            }
        }
        for (z = 0; z < dimz; z++) {
            exptr = &ex_psf[z*dimr];
            emptr = &em_psf[z*dimr];
            ovptr = &obsvol[z*dimr];
            /* emission psf in xy space */
            for (x = 0; x < dimd; x++) {
                indey = x;
                index = x*dimd;
                _em[index] = _em[indey] = emptr[x];
                for (y = 1; y <= x; y++) {
                    index++;
                    indey += dimd;
                    ri = _ri[index];
                    if (ri >= 0) {
                        rf = _rf[index];
                        _em[index] = _em[indey] =
                            rf ? emptr[ri]+rf*(emptr[ri+1]-emptr[ri])
                               : emptr[ri];
                    } else {
                        _em[index] = _em[indey] = 0.0;
                    }
                }
            }
            for (x = dimd; x < dimr; x++) {
                index = x*dimd;
                _em[index] = emptr[x];
                for (y = 1; y < dimd; y++) {
                    index++;
                    ri = _ri[index];
                    if (ri >= 0) {
                        rf = _rf[index];
                        _em[index] = rf ? emptr[ri]+rf*(emptr[ri+1]-emptr[ri])
                                        : emptr[ri];
                    } else {
                        _em[index] = 0.0;
                    }
                }
            }
            for (r = 0; r < dimr; r++) {
                /* Convolute emission PSF with detector kernel.
                For large kernels this is inefficient and should be
                replaced by a FFT based algorithm. */
                sum = 0.0;
                i = 1-dimd + (_dimd-dimd);
                for (x = r-dimd+1; x < std::min(r+dimd, dimr); x++) {
                    xx = abs(x) * dimd;
                    ii = abs(i++) * _dimd;
                    for (y = 0; y < dimd; y++) {
                        sum += _em[xx++] * detector[ii++];
                    }
                }
                /* multiply integral with excitation psf */
                ovptr[r] = sum * exptr[r];
            }
        }
        free(_ri);
        free(_rf);
        free(_em);
    }

    /* normalize maximum intensity */
    t = obsvol[0];
    for (i = 0; i < dimz*dimr; i++)
        obsvol[i] /= t;

    return 0;
}

/*
 *Calculate the detector kernel for integration over pinhole with trapezoid rule.
 *
 *The radius denotes the outer radius, except for square shape, where it denotes
 *the inner radius.
 */
int pinhole_kernel(int corners, double* out, int dim, double radius)
{
    int i, j, k;
    double alpha, t;

    for (i = 0; i < dim*dim; i++)
        out[i] = 1.0;

    switch (corners) {
        case 0: /* round pinhole */
            /* fill corner */
            t = sqrt(2.0 * dim*dim) - dim;
            t = t / sqrt(2.0);
            k = dim - ceil_int(t);
            for (i = k; i < dim; i++)
                for (j = k; j < dim; j++)
                    out[i*dim+j] = 0.0;
            /* draw antialiased arc using eightfold symmetry */
            for (j = 0; j <= floor_int((dim-1)/sqrt(2)); j++) {
                k = ceil_int(sqrt(radius*radius - (double)(j*j)));
                alpha = 1.0 - (sqrt((double)(k*k + j*j)) - radius);
                alpha *= 0.5;
                out[k*dim+j] = out[j*dim+k] = alpha;
                if (k > 0) {
                    k--;
                    alpha = radius - sqrt((double)(k*k + j*j));
                    alpha = 0.5 + 0.5 * alpha;
                    out[k*dim+j] = out[j*dim+k] = alpha;
                }
                for (i = k+2; i < dim; i++) {
                    out[i*dim+j] = out[j*dim+i] = 0.0;
                }
            }
            break;
        case 4: /* square pinhole */
            alpha = 0.5 * (radius - (double)floor_int(radius));
            for (i = 0; i < dim; i++) {
                out[i*dim + dim-1] *= alpha;
                out[i + (dim-1)*dim] *= alpha;
            }
            alpha = 0.5 + alpha;
            for (i = 0; i < dim-1; i++) {
                out[i*dim + dim-2] *= alpha;
                out[i + (dim-2)*dim] *= alpha;
            }
            break;
        default:
            return -1;
    }
    for (i = 0; i < dim; i++)
        for (j = 1; j < dim; j++)
            out[i*dim+j] *= 2.0;

    return 0;
}

/* Apply rotational symmetry around the 1st dimension.
 */
int zr2zxy(double* data, double* out, int dimz, int dimr)
{
    int x, y, z, x2, ri, index, indey, *_ri;
    double rf, *_rf, *dptr, *optr;

    /* floor integer and remainder float fraction of radius at xy */
    _ri = (int *)malloc((dimr*dimr)*sizeof(int));
    if (_ri == NULL)
        return -1;
    _rf = (double *)malloc((dimr*dimr)*sizeof(double));
    if (_rf == NULL) {
        free(_ri);
        return -1;
    }

    for (x = 0; x < dimr; x++) {
        x2 = x*x;
        indey = x;
        index = x*dimr;
        for (y = 0; y <= x; y++) {
            rf = sqrt((double)(x2 + y*y));
            ri = floor_int(rf);
            _ri[index] = _ri[indey] = (dimr > ri) ? ri : -1;
            _rf[index] = _rf[indey] = (dimr > ri+1) ? rf-(double)ri : 0.0;
            index++;
            indey += dimr;
        }
    }
    for (z = 0; z < dimz; z++) {
        dptr = &data[z*dimr];
        optr = &out[z*dimr*dimr];
        for (x = 0; x < dimr; x++) {
            indey = x;
            index = x*dimr;
            optr[index] = optr[indey] = dptr[x];
            for (y = 1; y <= x; y++) {
                index++;
                indey += dimr;
                ri = _ri[index];
                if (ri >= 0) {
                    rf = _rf[index];
                    optr[index] = optr[indey] =
                        rf ? dptr[ri]+rf*(dptr[ri+1]-dptr[ri]) : dptr[ri];
                } else {
                    optr[index] = optr[indey] = 0.0;
                }
            }
        }
    }
    free(_ri);
    free(_rf);
    return 0;
}

void psf_generate()
{
	
}
