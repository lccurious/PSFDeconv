#include <opencv2/opencv.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <Eigen/Eigen>
#include <cstdio>
#include <iostream>
#include "conv_img.hpp"

using namespace std;
using namespace cv;

Mat raw_image, out_img;
Mat kernel = Mat(3, 3, CV_8U);
Rect2d ROI_rect;

int main(int argc, char **argv)
{
	if (argc == 2)
	{
		raw_image = imread(argv[1]);
	}
	else
	{
		printf("Usage:\n\t%s %s\n", argv[0], "path_to_image");
	}
	
	ROI_rect = selectROI(raw_image, fromCenter=false);

	cout << ROI_rect.x <<endl;

	printf("image size is %d\t%d\n", raw_image.size().height, raw_image.size().width);

	out_img = deconv_img(raw_image, 1, 1, kernel);
	
	imshow("Cell", raw_image);
	imshow("Out", out_img);

	waitKey(-1);

	return 0;
}
