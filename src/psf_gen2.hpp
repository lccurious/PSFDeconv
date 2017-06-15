/* This header include the born & wolf generate algorithm
 * Author: hzn
 */
#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <complex>
#include <cmath>
#include <thread>
#include <atomic>

#include <boost/lambda/lambda.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <opencv2/opencv.hpp>

#include "vtkTools.hpp"
#include "cvTools.hpp"

typedef boost::multiprecision::cpp_dec_float_50 float_type;

typedef boost::math::policies::policy<
	boost::math::policies::domain_error<boost::math::policies::ignore_error>,
	boost::math::policies::overflow_error<boost::math::policies::ignore_error>,
	boost::math::policies::underflow_error<boost::math::policies::ignore_error>,
	boost::math::policies::denorm_error<boost::math::policies::ignore_error>,
	boost::math::policies::pole_error<boost::math::policies::ignore_error>,
	boost::math::policies::evaluation_error<boost::math::policies::ignore_error>
> ignore_all_policy;

int bessel_series()
{
	std::vector<double> bessel_table, X, bess_t;
	double bessel_tmp, bessel_sum = 0;
	double v = 0;
	
	// Always using try catch for unpredicted error;
	try
	{
		for (int i = 0; i < 2000; i++)
		{
			bessel_tmp = boost::math::sph_bessel(0, v);
			v += 0.01;
			bessel_sum += bessel_tmp * 0.01;
			X.push_back(v);
			bess_t.push_back(bessel_tmp);
			bessel_table.push_back(bessel_sum);
		}
		std::thread vtk_t(vtk_2Dplot, X, bessel_table);
		std::thread cv_t(cv_draw_pie, bessel_table);
		std::cout << "Viusalization finished!" << std::endl;
		vtk_t.join();
		cv_t.join();
	}
	catch( std::exception ex )
	{
		std::cout << "Thrown exception " << ex.what() << std::endl;
	}

	return 0;
}

double born_wolf_point(double k, double NA, double n_i, int x, int y, int z)
{
	double const_ratio = k * NA / n_i * std::sqrt(x*x + y*y);
	double const_opd = -0.5 * k * z * NA*NA/n_i/n_i;
	
	std::complex<double>bess_sum(0.0, 0.0);
	double bess_tmp, v = 0.0;
	int num_p = 500;
	double delta_v = 1.0 / 500;
	std::complex<double>opd(0.0, v);
	
	// Always using try for unpredicted error;
	try
	{
		for (int i = 0; i <= num_p; i++)
		{
			bess_tmp = boost::math::sph_bessel(0, v*const_ratio);
			opd.imag(v*v * const_opd);
			bess_sum += bess_tmp * std::exp(opd) * v;
			v += delta_v;
		}
		// std::cout << "integration:\t" << std::abs(bess_sum) << std::endl;
	}
	catch ( std::exception ex)
	{
		std::cout <<  "Thrown exception " << ex.what() << std::endl;
	}

	return std::abs(bess_sum);
}
