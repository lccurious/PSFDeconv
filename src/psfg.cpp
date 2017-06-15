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

int main(void)
{
	std::cout << "Here to show something ablout the end!" << std::endl;

	std::vector<double> X,Y;

	for (int i = 0; i < 1500; i++)
	{
		Y.push_back(born_wolf_point(M_2PI/ex_wavelen, NA, refr_index, i, i, 0));
		X.push_back(i);
	}
	vtk_2Dplot(X, Y);

  return 0;
}
