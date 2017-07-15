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
#endif //PSF_DECONVPSF_H
