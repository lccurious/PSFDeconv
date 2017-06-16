#include <thread>
#include <atomic>

#include "psf_gen2.hpp"
#include "vtkTools.hpp"

#define WIDTH_PSF 32
#define HEIGHT_PSF 32

double psf_data[WIDTH_PSF*HEIGHT_PSF];
int shape[2] = {WIDTH_PSF, HEIGHT_PSF};
double udim[2] = {4.0, 4.0};
double NA = 1.2;
double refr_index = 1.333;
double ex_wavelen = 488;
double em_wavelen = 520;
double pinhole_radius = 0.55;
#define M_2PI (6.283185307179586476925286766559)

int main(int argc, char** argv)
{
	std::cout << "Here to show something ablout the end!" << std::endl;
	/*
	std::vector<double> X,Y;
	cv::Mat img = Mat::zeros(800, 800, CV_8U);
	uchar *p1, *p2;
	for (int i = 0; i < 800; i++)
	{
		p1 = img.ptr<uchar>(799-i);
		for (int j = 0; j <= i; j++)
		{
			p1[j] = (int)born_wolf_point(M_2PI/ex_wavelen,NA, refr_index, i, j, 10);
			p2 = img.ptr<uchar>(799-j);
			p2[i] = p1[j];
		}
		std::cout<< i+1 << " / 800 finised!" << std::endl; 
		Y.push_back(p1[0]);
		X.push_back(i);
	}
	std::thread vtk_2d(vtk_2Dplot, X, Y);
	vtk_2d.join();
	cv::imshow("PSF 1/4", img);
	cv::waitKey(-1);
	*/
	int num_p = 512;
	std::vector<std::vector<double> >psf_matrix(num_p);
	born_wolf(0, psf_matrix,M_2PI/ex_wavelen, NA, refr_index, num_p);
	try
	{
		std::ofstream bess_file("Bessel.txt");
		for (int i = 0; i < num_p; i++)
		{
			for (int j = 0; j < num_p; j++)
			{
				bess_file << psf_matrix[i][j] << "\t";
			}
			bess_file << std::endl;
		}
		bess_file.close();
	}
	catch (std::exception ex)
	{
		std::cout << "Thrown excetion : " << ex.what() << std::endl;
	}
	show2DVec(psf_matrix);
  return 0;
}
