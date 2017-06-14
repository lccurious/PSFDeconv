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
	std::vector<double> bessel_table, X;
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
			std::cout << bessel_tmp << '\t' << bessel_sum << '\n';
			X.push_back(v);
			bessel_table.push_back(bessel_sum);
		}
		std::thread vtk_t(vtk_2Dplot, X, bessel_table);
	
		std::cout << "Viusalization finished!" << std::endl;
		vtk_t.join();
	}
	catch( std::exception ex )
	{
		std::cout << "Thrown exception " << ex.what() << std::endl;
	}

	return 0;
}

