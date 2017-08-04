//
// Created by peo on 17-7-13.
//

#ifndef PSF_DECONVPSF_H
#define PSF_DECONVPSF_H

#include <vector>
#include <cassert>
#include <iostream>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <fftw3.h>

#define OUTPUT_DEBUG std::cout<<"Passed here"<<std::endl;

/**
 * utility for testing.
 * @param frame
 * @param winName
 */
void norm_show(cv::InputArray frame, const char* winName);

/**
 * FFTW fourier transform interface.
 * @param squareSize
 * @param in_img
 * @param out_img
 */
void dftplane(int squareSize,
              cv::InputArray in_img,
              cv::OutputArray out_img);

/**
 * FFTW inv fourier transform interface.
 * @param squareSize
 * @param in_img
 * @param out_img
 */
void idftplane(int squareSize,
              cv::InputArray in_img,
              cv::OutputArray out_img);
/**
 * convolution between the raw image and the convolution core, use OpenCV DFT
 * @param A
 * @param B
 * @param C
 */
void convolutionDFT(cv::InputArray A, cv::InputArray B, cv::OutputArray C);


/**
 * The using openCV dft to transform img into fourier presentation, and show the normalized
 * fourier domain image with openCV highgui
 * @param img the input image
 */
void fourier_show(cv::InputArray img);

/**
 * execute the FFT by using method in FFTW, Only process the one channel
 * @param in_img
 * @param out_img
 * @return
 */
int single_fft(cv::InputArray in_img, cv::OutputArray out_img);

/**
 * TODO(peo):Once CUDA got, change this function, before that, please do not use this function.
 * convolution operation between the Input image A and B
 * @param A
 * @param B
 * @param C
 */
void multiply_fourier(cv::InputArray A, cv::InputArray B, cv::OutputArray C);

/**
 * only accept the CV_64F type, respectively the double type
 * @param A Input image for PSF core for convolution in general
 * @param B PSF core in general
 * @param C Output image
 */
void divide_fourier(cv::InputArray A, cv::InputArray B, cv::OutputArray C);

/**
 * for testing frequency matrix FFT transform
 * @param A
 * @param B
 * @param C
 */
void divide_fourier_test(cv::InputArray A, cv::InputArray B, cv::OutputArray C);

/**
 * The shape of _srcA and _srcB must have the same shape and data type
 * @param _srcA
 * @param _srcB
 * @param _dst
 * @param flags
 * @param conjB
 */
static void divSpectrums(cv::InputArray _srcA, cv::InputArray _srcB, cv::OutputArray _dst, int flags, bool conjB);

/**
 * Richardson-Lucy implementation
 * @param _srcI
 * @param _coreI
 * @param _dst
 * @param iteration
 */
void RichardLucy(cv::InputArray _srcI,
                 cv::InputArray _coreI,
                 cv::OutputArray _dst,
                 int iteration);

/**
 * testing the way padded method influence the RL-method.
 * @param _srcI
 * @param _coreI
 * @param _dst
 */
void RichardLucy_single(cv::InputArray _srcI, cv::InputArray _coreI, cv::OutputArray _dst);

#endif //PSF_DECONVPSF_H
