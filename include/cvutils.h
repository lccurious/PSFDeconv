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

#define INF 65535

int cv_draw_pie(std::vector<double>amp);
int simple_show(const std::vector<std::vector<double> >&plane);
int show2DVec(const std::vector<std::vector<double> >&plane);
int getTIFF(const char* filename, cv::OutputArray);
void vec2mat(std::vector<std::vector<double> > M2D, const cv::OutputArray out_im);

#endif //PSF_CVUTILS_H
