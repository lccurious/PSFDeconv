#include <opencv2/opencv.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <Eigen/Eigen>
#include <cstdio>
#include <iostream>
#include "conv_img.hpp"
#include "wiener_deconv.hpp"

using namespace std;
using namespace cv;

Mat raw_image, out_img;
Rect2d ROI_rect;

// Define the the kernel.
Mat kernel = Mat::ones(3, 3, CV_32F) / (float)(9);

int main(int argc, char **argv)
{
	if (argc == 2)
	{
		raw_image = imread(argv[1], 0);
	}
	else
	{
		printf("Usage:\n\t%s %s\n", argv[0], "path_to_image");
	}
	
	printf("image size is %d\t%d\n", raw_image.size().height, raw_image.size().width);
	
	imshow("Cell", raw_image);
	imshow("Out", out_img);

	waitKey(-1);

	return 0;
}
