/* This header include the born & wolf generate algorithm
 * Author: hzn
 */

#ifndef PSFGEN_H
#define PSFGEN_H

#include <vector>
#include <iterator>
#include <complex>
#include <cmath>

#include <boost/lambda/lambda.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/progress.hpp>
#include <opencv2/opencv.hpp>

extern int AiryRadius;

typedef boost::multiprecision::cpp_dec_float_50 float_type;

typedef boost::math::policies::policy<
	boost::math::policies::domain_error<boost::math::policies::ignore_error>,
	boost::math::policies::overflow_error<boost::math::policies::ignore_error>,
	boost::math::policies::underflow_error<boost::math::policies::ignore_error>,
	boost::math::policies::denorm_error<boost::math::policies::ignore_error>,
	boost::math::policies::pole_error<boost::math::policies::ignore_error>,
	boost::math::policies::evaluation_error<boost::math::policies::ignore_error>
> ignore_all_policy;

/* 一个测试函数，用于测试贝塞尔函数是否能够正常生成
 */
int bessel_series(void);

/* 对于某一个单独的像素点，生成基于born&wolf点扩散函数在那个像素的值
 */
double born_wolf_point(double k, double NA, double n_i, int x, int y, int z);

/* 对于某一个z值，生成在对应的整个平面的点扩散函数二维矩阵
 */
int born_wolf(int z, std::vector<std::vector<double> >& M2D,
              double k, double NA, double n_i, int num_p);

/**
 * build an complete PSF plane for a given z
 * @param z
 * @param M2D
 * @param k
 * @param NA
 * @param n_i
 * @param num_p
 * @return
 */
int born_wolf_full(int z, std::vector<std::vector<double> >& M2D,
				   double k, double NA, double n_i, int num_p);


/* Test function for bessel integral
 */
int born_wolf_test(int z, std::vector<std::vector<double> >& M2D,
				   double k, double NA, double n_i, int num_p);

/**!
 * only for test the integral operation.
 * @param z the value in z axial
 * @param ori_vec original vector
 * @param int_vec
 * @param num_p
 * @param xSize
 * @return
 */
int integral_test(int z,
                  std::vector<double>& ori_vec,
                  std::vector<double>& int_vec,
                  int num_p,
                  int xSize);

/**
 * with full stack generation;
 * @param PSF3D
 * @param numstack
 * @param focal_idx
 * @param PSFSize
 * @param k
 * @param NA
 * @param n_i
 * @return
 */
int BornWolf_stack(std::vector<cv::Mat> &PSF3D, int numstack, int focal_idx,
                   int PSFSize,
                   double PIK, double NA, double n_i);

#endif //PSFGENH