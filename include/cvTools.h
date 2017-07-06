//
// Created by peo on 17-6-27.
//

#ifndef PSF_CVTOOLS_H
#define PSF_CVTOOLS_H

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

#define INF 65535

/* 对于当前的
 */
int cv_draw_pie(std::vector<double>amp);

int show2DVec(const std::vector<std::vector<double> >&plane);

#endif //PSF_CVTOOLS_H
