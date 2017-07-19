//
// Created by peo on 17-7-10.
//

#ifndef PSF_CVUTILS_H
#define PSF_CVUTILS_H

#include <iostream>
#include <fstream>
#include <tiffio.h>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#define INF 65535

int cv_draw_pie(std::vector<double>amp);
// round the 1D vector into circular plane

int simple_show_vec(const std::vector<std::vector<double> >&plane);
// show plane for observation only

int show2DVec(const std::vector<std::vector<double> >&plane);
// visualization of single 2D vector

int getTIFFinfo(const char* filename);
// print some TIFF Image information

int TIFFframenumber(const char* filename);
// count the entire TIFF Image frame number

int getTIFF(const char* filename, cv::OutputArray out_img, uint16 frame_seq);
// out_im with CV_64F data type

void vec2mat(std::vector<std::vector<double> > M2D, const cv::OutputArray out_im);
// out_im with CV_64F data type

int reconstructZY(const char* filename, cv::OutputArray out_img);

static void on_trackbar(int , void*);

int cv_slideWin(const char* filename);

#endif //PSF_CVUTILS_H
