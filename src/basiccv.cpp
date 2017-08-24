//
// Created by peo on 17-8-10.
//

#include "basiccv.h"


void high_reject_filter(cv::InputArray _srcI)
{
    cv::Mat srcI = _srcI.getMat();
    cv::Mat padded;
    int opt_width = cv::getOptimalDFTSize(srcI.rows);
    int opt_height = cv::getOptimalDFTSize(srcI.cols);
    cv::copyMakeBorder(srcI, padded, 0, opt_width-srcI.rows, 0, opt_height-srcI.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));

    // DFT
    cv::Mat planes[] = {cv::Mat_<double>(padded), cv::Mat::zeros(padded.size(), CV_64F)};
    cv::Mat com_img;
    cv::merge(planes, 2, com_img);
    cv::dft(com_img, com_img);

    // get DFT image
    cv::split(com_img, planes);
    cv::magnitude(planes[0], planes[1], planes[0]);

    cv::Mat mag_mat = planes[0];
    mag_mat += cv::Scalar::all(1);
    cv::log(mag_mat, mag_mat);

    // relocation quad image
    mag_mat = mag_mat(cv::Rect(0, 0, mag_mat.cols & -2, mag_mat.rows & -2));
    int cx = mag_mat.cols / 2;
    int cy = mag_mat.rows / 2;

    cv::Mat q0(mag_mat, cv::Rect(0, 0, cx, cy));
    cv::Mat q1(mag_mat, cv::Rect(0, cy, cx, cy));
    cv::Mat q2(mag_mat, cv::Rect(cx, cy, cx, cy));
    cv::Mat q3(mag_mat, cv::Rect(cx, 0, cx, cy));

    cv::Mat tmp;
    q0.copyTo(tmp);
    q2.copyTo(q0);
    tmp.copyTo(q2);

    q1.copyTo(tmp);
    q3.copyTo(q1);
    tmp.copyTo(q3);



    // prepare for exhibition
    cv::Mat mag_img(mag_mat);
}