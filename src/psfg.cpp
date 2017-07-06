#include "cvTools.h"
#include "vtkTools.h"
#include "psf_gen2.h"

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
	clock_t start, finish;
	double duration;
	std::cout << "Here to show something ablout the end!" << std::endl;
	int num_p = 256;
    int stack_depth = 1;

	// 预分配程序矩阵空间
	std::vector<std::vector<double> >psf_matrix(num_p);
	std::vector<std::vector<double> >psf_matrix_copy(num_p);
	boost::progress_display *show_progress = NULL;
    show_progress = new boost::progress_display(stack_depth);
	
	// 性能测试计时开始
    start = clock();
    for (int i = 0; i < stack_depth; i++)
	{
		born_wolf(i, psf_matrix, M_2PI/ex_wavelen, NA, refr_index, num_p);
        ++(*show_progress);
	}
	
	// 性能测试计时结束
    finish = clock();
    std::cout << "Total Time spent: " << (finish - start)/CLOCKS_PER_SEC
              << " Seconds" << std::endl;

	show2DVec(psf_matrix);
    return 0;
}
