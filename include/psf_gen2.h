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

int bessel_series(void);

double born_wolf_point(double k, double NA, double n_i, int x, int y, int z);

int born_wolf(int z, std::vector<std::vector<double> >& M2D,
              double k, double NA, double n_i, int num_p);

int born_wolf_zstack(std::vector<int> zs,
                     std::vector<std::vector<std::vector<double> > >& M3D,
                     double k, double NA, double n_i, int num_p);

#endif //PSF_GEN2_H