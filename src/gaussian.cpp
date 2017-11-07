/**
 * \file gaussian_sii.cpp
 * \brief Gaussian convolution using stacked integral images
 * \author Pascal Getreuer <getreuer@cmla.ens-cachan.fr>
 *
 * Copyright (c) 2012-2013, Pascal Getreuer
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under, at your option, the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version, or the terms of the
 * simplified BSD license.
 *
 * You should have received a copy of these licenses along with this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include "gaussian.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
/** \brief The constant pi */
#define M_PI 3.14159265358979323846264338327950288
#endif


/** \brief Minimum SII filter order */
#define SII_MIN_K       3
/** \brief Maximum SII filter order */
#define SII_MAX_K       5
/** \brief Test whether a given K value is a valid SII filter order */
#define SII_VALID_K(K)  (SII_MIN_K <= (K) && (K) <= SII_MAX_K)


/** \brief Parameters for stacked integral images Gaussian approximation */
typedef struct sii_coeffs_
{
    float weights[SII_MAX_K];   /**< Box weights     */
    long radii[SII_MAX_K];      /**< Box radii       */
    int K;                      /**< Number of boxes */
} sii_coeffs;


/**
 * \brief Precompute filter coefficients for SII Gaussian convolution
 * \param c         sii_coeffs pointer to hold precomputed coefficients
 * \param sigma     Gaussian standard deviation
 * \param K         number of boxes = 3, 4, or 5
 * \return 1 on success, 0 on failure
 *
 * This routine reads Elboher and Werman's optimal SII radii and weights for
 * reference standard deviation \f$ \sigma_0 = 100/\pi \f$ from a table and
 * scales them to the specified value of sigma.
 */
void sii_precomp(sii_coeffs *c, double sigma, int K)
{
    /* Elboher and Werman's optimal radii and weights. */
    const double sigma0 = 100.0 / M_PI;
    static const short radii0[SII_MAX_K - SII_MIN_K + 1][SII_MAX_K] =
        {{76, 46, 23, 0, 0},
         {82, 56, 37, 19, 0},
         {85, 61, 44, 30, 16}};
    static const float weights0[SII_MAX_K - SII_MIN_K + 1][SII_MAX_K] =
        {{0.1618f, 0.5502f, 0.9495f, 0, 0},
         {0.0976f, 0.3376f, 0.6700f, 0.9649f, 0},
         {0.0739f, 0.2534f, 0.5031f, 0.7596f, 0.9738f}};

    const int i = K - SII_MIN_K;
    double sum;
    int k;
    
    //assert(c && sigma > 0 && SII_VALID_K(K));
    c->K = K;
    
    for (k = 0, sum = 0; k < K; ++k)
    {
        c->radii[k] = (long)(radii0[i][k] * (sigma / sigma0) + 0.5);
        sum += weights0[i][k] * (2 * c->radii[k] + 1);
    }
    
    for (k = 0; k < K; ++k)
        c->weights[k] = (float)(weights0[i][k] / sum);
    
    return;
}

/**
 * \brief Determines the buffer size needed for SII Gaussian convolution
 * \param c     sii_coeffs created by sii_precomp()
 * \param N     number of samples
 * \return required buffer size in units of num samples
 *
 * This routine determines the minimum size of the buffer needed for use in
 * sii_gaussian_conv() or sii_gaussian_conv_image(). This size is the length
 * of the signal (or in 2D, max(width, height)) plus the twice largest box
 * radius, for padding.
 */
long sii_buffer_size(sii_coeffs c, long N)
{
    long pad = c.radii[0] + 1;
    return N + 2 * pad;
}



/**
 * \brief Half-sample symmetric boundary extension
 * \param N     signal length
 * \param n     requested sample, possibly outside {0,...,`N`-1}
 * \return reflected sample in {0,...,`N`-1}
 *
 * This function is used for boundary handling. Suppose that `src` is an array
 * of length `N`, then `src[extension(N, n)]` evaluates the symmetric
 * extension of `src` at location `n`.
 *
 * Half-sample symmetric extension is implemented by the pseudocode
\verbatim
    repeat
        if n < 0, reflect n over -1/2
        if n >= N, reflect n over N - 1/2
    until 0 <= n < N
\endverbatim
 * The loop is necessary as some `n` require multiple reflections to bring
 * them into the domain {0,...,`N`-1}.
 *
 * This function is used by all of the Gaussian convolution algorithms
 * included in this work except for DCT-based convolution (where symmetric
 * boundary handling is performed implicitly by the transform). For FIR, box,
 * extended box, SII, and Deriche filtering, this function could be replaced
 * to apply some other boundary extension (e.g., periodic or constant
 * extrapolation) without any further changes. However, Alvarez-Mazorra and
 * Vliet-Young-Verbeek are hard-coded for symmetric extension on the right
 * boundary, and would require specific modification to change the handling
 * on the right boundary.
 *
 * \par A note on efficiency
 * This function is a computational bottleneck, as it is used within core
 * filtering loops. As a small optimization, we encourage inlining by defining
 * the function as `static`. We refrain from further optimization since this
 * is a pedagogical implementation, and code readability is more important.
 * Ideally, filtering routines should take advantage of algorithm-specific
 * properties such as exploiting sequential sample locations (to update the
 * extension cheaply) and samples that are provably in the interior (where
 * boundary checks may omitted be entirely).
 */
static long extension(long N, long n)
{
    while (1)
        if (n < 0)
            n = -1 - n;         /* Reflect over n = -1/2.    */
        else if (n >= N)
            n = 2 * N - 1 - n;  /* Reflect over n = N - 1/2. */
        else
            break;
        
    return n;
}



/**
 * \brief Gaussian convolution SII approximation
 * \param c         sii_coeffs created by sii_precomp()
 * \param dest      output convolved data
 * \param buffer    array with space for sii_buffer_size() samples
 * \param src       input, modified in-place if src = dest
 * \param N         number of samples
 * \param stride    stride between successive samples
 *
 * This routine performs stacked integral images approximation of Gaussian
 * convolution with half-sample symmetric boundary handling. The buffer array
 * is used to store the cumulative sum of the input, and must have space for
 * at least sii_buffer_size(c,N) samples.
 *
 * The convolution can be performed in-place by setting `src` = `dest` (the
 * source array is overwritten with the result). However, the buffer array
 * must be distinct from `src` and `dest`.
 */
void sii_gaussian_conv(sii_coeffs c, float *dest, float *buffer,
    const float *src, long N, long stride)
{
    float accum;
    long pad, n;
    int k;
    
    //assert(dest && buffer && src && dest != buffer && src != buffer
    //    && N > 0 && stride != 0);
    
    pad = c.radii[0] + 1;
    buffer += pad;
    
    /* Compute cumulative sum of src over n = -pad,..., N + pad - 1. */
    for (n = -pad, accum = 0; n < N + pad; ++n)
    {
        accum += src[stride * extension(N, n)];
        buffer[n] = accum;
    }
    
    /* Compute stacked box filters. */
    for (n = 0; n < N; ++n, dest += stride)
    {
        accum = c.weights[0] * (buffer[n + c.radii[0]]
            - buffer[n - c.radii[0] - 1]);
        
        for (k = 1; k < c.K; ++k)
            accum += c.weights[k] * (buffer[n + c.radii[k]]
                - buffer[n - c.radii[k] - 1]);
        
        *dest = accum;
    }
    
    return;
}

/**
 * \brief 2D Gaussian convolution SII approximation
 * \param c             sii_coeffs created by sii_precomp()
 * \param dest          output convolved data
 * \param buffer        array with space for sii_buffer_size() samples
 * \param src           image to be convolved, overwritten if src = dest
 * \param width         image width
 * \param height        image height
 * \param num_channels  number of image channels
 *
 * Similar to sii_gaussian_conv(), this routine approximates 2D Gaussian
 * convolution with stacked integral images. The buffer array must have space
 * for at least sii_buffer_size(c,max(width,height)) samples.
 *
 * The convolution can be performed in-place by setting `src` = `dest` (the
 * source array is overwritten with the result). However, the buffer array
 * must be distinct from `src` and `dest`.
 */
void sii_gaussian_conv_image(sii_coeffs c, float *dest, float *buffer,
    const float *src, int width, int height, int num_channels)
{
    long num_pixels = ((long)width) * ((long)height);
    int x, y, channel;
    
    //assert(dest && buffer && src && num_pixels > 0);
    
    /* Loop over the image channels. */
    for (channel = 0; channel < num_channels; ++channel)
    {
        float *dest_y = dest;
        const float *src_y = src;
        
        /* Filter each row of the channel. */
        for (y = 0; y < height; ++y)
        {
            sii_gaussian_conv(c,
                dest_y, buffer, src_y, width, 1);
            dest_y += width;
            src_y += width;
        }
        
        /* Filter each column of the channel. */
        for (x = 0; x < width; ++x)
            sii_gaussian_conv(c,
                dest + x, buffer, dest + x, height, width);
        
        dest += num_pixels;
        src += num_pixels;
    }
    
    return;
}


/**
 *
 * Convolution with a Gaussian using Stacked Integral Images
 *
 */
void gaussian(
  float *I,     //input/output image
  int   nx,     //image width
  int   ny,     //image height
  float sigma,  //Gaussian sigma
  int   K       //defines the number of iterations
)
{
  float *buffer = NULL;
  sii_coeffs c;
    
  sii_precomp(&c, sigma, K);
  
  if (!(buffer = (float *)malloc(sizeof(float) * sii_buffer_size(c,
      ((nx >= ny) ? nx : ny)))))
      return;
  
  sii_gaussian_conv_image(c, I, buffer, I, nx, ny, 1);
  free(buffer);
}

