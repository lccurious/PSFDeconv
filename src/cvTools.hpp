#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

#define INF 65535

int cv_draw_pie(std::vector<double>amp)
{
	int width, height, cx, cy;
	int num_v = amp.size();
	width = height = amp.size() * 2 + 10;
	cx = width / 2;
	cy = height / 2;
	cv::Mat img;
	img.create(width, height, 1);
	
	double amp_min = INF;
	double amp_max = -INF;

	for (std::vector<double>::iterator iter = amp.begin(); iter < amp.end(); iter++)
	{
		if (*iter < amp_min)
			amp_min = *iter;
		if (amp_min < *iter)
			amp_max = *iter;
	}
	double differ_dis = amp_max - amp_min;

	for (int i = 0; i < num_v; i++)
	{
		circle(img, cv::Point(cx, cy), i, (int)((amp[i]-amp_min)*255 / differ_dis));
	}

	if (img.rows > 600)
	{
		cv::Mat dst;
		double fx = (double)600 / img.cols;
		resize(img, dst,cv::Size(), fx, fx);
		imshow("PSF Plane", dst);
	}
	else
	{
		imshow("PSF Plane", img);
	}
	
	cv::waitKey(-1);

	return 0;
}


int show2DVec(const std::vector<std::vector<double> >&plane)
{
	if (plane.size() < 1)
	{
		std::cout << "Plane size not leggal" << std::endl;
	}
	int width, height;
	height = plane.size();
	width = height;
	std::cout << "I got the Width: " << width << "\tHeight: " << height << std::endl;
	cv::Mat img = cv::Mat::zeros(height*2, width*2, CV_8U);
	uchar *p1, *p2;
	double max_pixel = plane[0][0];
	for (int i = 0; i < height; i++)
	{
		p1 = img.ptr<uchar>(height + i);
		p2 = img.ptr<uchar>(height - i);
		for (int j = 0; j < width; j++)
		{
			p1[width+j] = plane[i][j] * 255 / max_pixel;
			p1[width-j] = plane[i][j] * 255 / max_pixel;
			p2[width+j] = plane[i][j] * 255 / max_pixel;
			p2[width-j] = plane[i][j] * 255 / max_pixel;
		}
	}
	cv::imshow("Complete PSF", img);
	cv::imwrite("images/Complete PSF.png", img);
	cv::waitKey(-1);
	return 0;
}
