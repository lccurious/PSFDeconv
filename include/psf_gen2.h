/* This header include the born & wolf generate algorithm
 * Author: hzn
 */

#ifndef PSF_GEN2_H
#define PSF_GEN2_H

#include <vector>
#include <iterator>
#include <complex>
#include <cmath>

#include <boost/lambda/lambda.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/progress.hpp>


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

/* 生成一个完整的点扩散函数的矩阵
 */
int born_wolf_full(int z, std::vector<std::vector<double> >& M2D,
				   double k, double NA, double n_i, int num_p);

int born_wolf_zstack(std::vector<int> zs,
                     std::vector<std::vector<std::vector<double> > >& M3D,
                     double k, double NA, double n_i, int num_p);

/* Test function for bessel integral
 */
int born_wolf_test(int z, std::vector<std::vector<double> >& M2D,
				   double k, double NA, double n_i, int num_p);

int integral_test(int z,
                  std::vector<double>& ori_vec,
                  std::vector<double>& int_vec,
                  int num_p,
                  int xSize);

#endif //PSF_GEN2_H