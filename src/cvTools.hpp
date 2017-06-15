#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

using namespace std;
using namespace cv;

#define INF 65535

int cv_draw_pie(vector<double>amp)
{
	int width, height, cx, cy;
	int num_v = amp.size();
	width = height = amp.size() * 2 + 10;
	cx = width / 2;
	cy = height / 2;
	Mat img;
	img.create(width, height, 1);
	
	double amp_min = INF;
	double amp_max = -INF;

	for (vector<double>::iterator iter = amp.begin(); iter < amp.end(); iter++)
	{
		if (*iter < amp_min)
			amp_min = *iter;
		if (amp_min < *iter)
			amp_max = *iter;
	}
	double differ_dis = amp_max - amp_min;

	for (int i = 0; i < num_v; i++)
	{
		circle(img, Point(cx, cy), i, (int)((amp[i]-amp_min)*255 / differ_dis));
	}

	if (img.rows > 600)
	{
		Mat dst;
		double fx = (double)600 / img.cols;
		resize(img, dst, Size(), fx, fx);
		imshow("PSF Plane", dst);
	}
	else
	{
		imshow("PSF Plane", img);
	}
	
	waitKey(-1);

	return 0;
}
