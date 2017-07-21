//
// Created by peo on 17-7-13.
//

#ifndef PSF_DECONVPSF_H
#define PSF_DECONVPSF_H
#include <vector>
#include <cassert>
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/core/eigen.hpp>
#include <fftw3.h>

#define OUTPUT_DEBUG std::cout<<"Passed here"<<std::endl;

/**
 * I forget it.
 * @param raw_im
 * @param dis_im
 * @return
 */
int deconv_im(const std::vector<std::vector<double> > raw_im,
              std::vector<std::vector<double> >&dis_im);

/**
 * I forget it.
 * @param in_im
 * @param out_im
 * @return
 */
int cv_deconv_im(cv::InputArray in_im, cv::OutputArray out_im);

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
void iftplane(int squareSize,
              cv::InputArray in_img,
              cv::OutputArray out_img);
/**
 * convolution between the raw image and the convolution core
 * @param A
 * @param B
 * @param C
 */
void convolutionDFT(cv::InputArray A, cv::InputArray B, cv::OutputArray C);

/**
 * Richardson-Lucy reconstruction
 * @param img input image
 * @param core construct PSF core
 * @param out_img the restoration result
 */
void RichardLucydeconv(cv::InputArray img,
                       cv::InputArray core,
                       cv::OutputArray out_img);

/**
 * The using openCV dft to transform img into fourier presentation, and show the normalized
 * fourier domain image with openCV highgui
 * @param img the input image
 */
void fourier_show(cv::InputArray img);

int single_fft(cv::InputArray in_img, cv::OutputArray out_img);

// TODO(peo):Once CUDA got, change this function, before that, please do not use this function.
void multiply_fourier(cv::InputArray A, cv::InputArray B, cv::OutputArray C);

/**
 * only accept the CV_64F type, respectively the double type
 * @param A
 * @param B
 * @param C
 */
void divide_fourier(cv::InputArray A, cv::InputArray B, cv::OutputArray C);

// The shape of _srcA and _srcB must be the same shape and data type
static void divSpectrums(cv::InputArray _srcA, cv::InputArray _srcB, cv::OutputArray _dst, int flags, bool conjB);

#endif //PSF_DECONVPSF_H
