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

int deconv_im(const std::vector<std::vector<double> > raw_im,
              std::vector<std::vector<double> >&dis_im);
int cv_deconv_im(cv::InputArray in_im, cv::OutputArray out_im);

void dftplane(int squareSize,
              cv::InputArray in_img,
              cv::OutputArray out_img);

void iftplane(int squareSize,
              cv::InputArray in_img,
              cv::OutputArray out_img);

void convolutionDFT(cv::InputArray A, cv::InputArray B, cv::OutputArray C);

void RichardLucydeconv(cv::InputArray img,
                       cv::InputArray core,
                       cv::OutputArray out_img);

void fourier_show(cv::InputArray img);

int single_fft(cv::InputArray in_img, cv::OutputArray out_img);

void multiply_fourier(cv::InputArray A, cv::InputArray B, cv::OutputArray C);

static void divSpectrums(cv::InputArray _srcA, cv::InputArray _srcB, cv::OutputArray _dst, int flags, bool conjB);

#endif //PSF_DECONVPSF_H
