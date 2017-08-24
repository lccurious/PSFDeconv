//
// Created by peo on 17-8-10.
//

#ifndef PSF_BASICCV_H
#define PSF_BASICCV_H
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

/**
 * High frequency reject filter. try to filter the highfrequency noise in image.
 * @param _srcI Input Image, will be subsitatuted
 */
void high_reject_filter(cv::InputArray _srcI);

#endif //PSF_BASICCV_H
