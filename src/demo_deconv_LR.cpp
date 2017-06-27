//
// Created by peo on 17-6-27.
//

#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

#define kappa 10000

int main()
{
    cv::Mat im;
    cv::Mat im_conv_kernel;
    cv::Mat im_correction;
    cv::Mat im_new;
    cv::Mat im_new_est;
    cv::Mat im1;

    char str[80];
    im1 = cv::imread("images/B0008393-800.jpg");
    cv::imshow("Original-Color", im1);

    im = cv::imread("images/B0008393-800.jpg", 0);
    cv::imshow("Original-Gray", im);

    for (int i = 0; i < 100; i++) {

        double a[9] = {0, 40, 0, 0, 40, 0, 0, 40, 0};
        cv::Mat kernel = cv::Mat(3, 3, CV_32FC1, a);

        cv::filter2D(im, im_conv_kernel, im.depth(), kernel, cv::Point(-1, -1));

        cv::imshow("conv", im_conv_kernel);

        // Subtract from blurred image. Error correction = b(x, y) - ik(x, y)**k(x, y)
        im.copyTo(im_correction);
        cv::subtract(im, im_conv_kernel, im_correction);
        cv::imshow("Sub", im_correction);

        // Add ik(x, y) with image Correction - ik(x,y) + b(x,y) - ik(x,y) ** k(x,y)
        im.copyTo(im_new_est);
        cv::add(im, im_correction, im_new_est);
        cv::imshow("Add", im_new_est);

        im_new_est.copyTo(im);
    }
    cvWaitKey(-1);
    return 0;
}