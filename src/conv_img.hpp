#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/opencv.hpp>

using namespace cv;
using namespace std;

Mat deconv_img(Mat src, int padding, int strides, Mat kernel)
{

	int channels = src.channels();
	int nRows = src.rows;
	int nCols = src.cols * channels;
	int kernel_size = kernel.rows;

	int im_h = src.rows + padding*2;
	int im_w = src.cols + padding*2;
	im_w = im_w - kernel.cols / 2;
	im_h = im_h - kernel.rows / 2;
	im_w /= strides;
	im_h /= strides;
	printf("The im_h: %d\tim_w: %d\tkernel: %d\n", im_h, im_w, kernel.rows);
	
	Mat out_im(im_h, im_w, CV_8U);
	uchar *p;
	uchar *im_p;
	for (int i = 0; i < im_h; i++)
	{
		im_p = out_im.ptr<uchar>(i);
		p = src.ptr<uchar>(i);
		for (int j = 0; j < nCols; j++)
		{
			im_p[j] = p[j*3];
			// for (int k1 = 0; k1 < kernel_size; k1++)
			// {
			// 	p = src.ptr<uchar>(i);
			// 	for (int k2 = 0; k2 < kernel_size; k2++)
			// 	{
			// 		im_p[j] = p[k2];
			// 	}
			// }
		}
	}
	return out_im; 
}
