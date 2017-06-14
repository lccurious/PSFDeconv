#include "psf_gen2.hpp"

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

int main(void)
{
	std::cout << "Here to show something ablout the end!" << std::endl;

  bessel_series();

  return 0;
}
