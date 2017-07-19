//
// Created by peo on 17-6-27.
//

/* This header include the born & wolf generate algorithm
 * Author: hzn
 */
#include "psf_gen2.h"

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
    // 这是用于控制积分精度的，从0积到1这个数值越大，则积分精度越高
    // 如果要追求绝对的精度的话，要做一个trade off，就是数值精度的问题
    int num_p = 1000;
    double delta_v = 1.0 / num_p;
    std::complex<double>opd(0.0, v);

    // Always using try for unpredicted error;
    try
    {
        for (int i = 0; i <= num_p; i++)
        {
            bess_tmp = boost::math::sph_bessel(0, v*const_ratio);

            // TODO(peo):这里还会有大改动，可能计算过程出现失误
            opd.imag(v*v * const_opd);
            bess_sum += bess_tmp * std::exp(opd) * delta_v;
            v += delta_v;
        }
    }
    catch ( std::exception ex)
    {
        std::cout <<  "Thrown exception " << ex.what() << std::endl;
    }

    return std::fabs(bess_sum);
}

int born_wolf(int z, std::vector<std::vector<double> >& M2D,
              double k, double NA, double n_i, int num_p)
{
    int step = 7000 / num_p;
    double bessel_res = 0.0;
    try
    {
        for (int i = 0; i < num_p; i++)
        {
            M2D[i].resize(num_p);
            for (int j = 0; j <= i; j++)
            {
                bessel_res = born_wolf_point(k, NA, n_i, j*step, i*step, z);
                M2D[i][j] = bessel_res;
                M2D[j][i] = bessel_res;
            }
        }
    }
    catch (std::exception ex)
    {
        std::cout << "Thrown exception : " << ex.what() << std::endl;
    }
    return 0;
}

int born_wolf_full(int z, std::vector<std::vector<double> >& M2D,
                   double k, double NA, double n_i, int num_p)
{
    std::vector<std::vector<double> >M2D_cp(num_p);
    double step = 7000 / num_p;
    double bessel_res = 0.0;
    M2D.resize(num_p*2);
#ifdef _OBSERVE_MAX_PIXEL
    std::cout << "Max_pixel : " << max_pixel << std::endl;
#endif
    try {
        for (int i = 0; i < num_p; i++) {
            M2D[i].resize(num_p*2);
            M2D[i+num_p].resize(num_p*2);
            M2D_cp[i].resize(num_p);
            for (int j = 0; j <= i; j++) {
                bessel_res = born_wolf_point(k, NA, n_i, j*step, i*step, z);
                M2D_cp[i][j] = bessel_res;
                M2D_cp[j][i] = bessel_res;
            }
        }
        for (int i = 0; i < num_p; ++i) {
            for (int j = 0; j < num_p; ++j) {
                M2D[i+num_p][j+num_p] = M2D_cp[i][j];
                M2D[num_p-i-1][j+num_p] = M2D_cp[i][j];
                M2D[i+num_p][num_p-j-1] = M2D_cp[i][j];
                M2D[num_p-i-1][num_p-j-1] = M2D_cp[i][j];
            }
        }
    }
    catch (std::exception ex)
    {
        std::cout << "Thrown exception : " << ex.what() << std::endl;
    }
    return 0;
}

// TODO(peo): Here complete the PSF accelerate algorithm.
//double born_wolf3D(double k, double NA, double refr_n, int psfSize, int stackSize)
//{
//    std::vector<std::vector<double> > baseM[psfSize];
//    double const_ratio, const_opd;
//    for (int x = 0; x < psfSize; ++x) {
//        baseM[x].resize(psfSize);
//        for (int y = 0; y < psfSize; ++y) {
//            const_ratio = k*NA / refr_n*sqrt(x*x+y*y);
//
//        }
//    }
//    double const_ratio = k * NA / refr_n * sqrt(x*x + y*y);
//}

int born_wolf_test(int z, std::vector<std::vector<double> >& M2D,
                   double k, double NA, double n_i, int num_p)
{
    std::complex<double> test_a;
    std::complex<double> test_res;
    double delta_v = 1.0 / num_p;
    double v = 0.0;
    try {
        for (int i = 0; i <= num_p; ++i) {
            test_a.imag(v*2);
            test_res = boost::math::sph_bessel(0, v) * exp(test_a);
            std::cout <<"exp(J_0("<<v<<")) = "
                      <<test_res
                      <<"abs(test_res):\t"
                      <<fabs(test_res)
                      <<std::endl;
            v += delta_v;
        }
    }
    catch (std::exception ex)
    {
        std::cout << "Thrown exception : " << ex.what() << std::endl;
    }

    return 0;
}

int integral_test(int z,
                  std::vector<double>& ori_vec,
                  std::vector<double>& int_vec,
                  int num_p,
                  int xSize)
{
    std::complex<double> test_ori, test_res;
    double delta_v = 10.0 / num_p;
    double v = 0.0;
    double bess_sum = 0., bess_tmp;
    try {
        test_ori.imag(0.);
        test_ori.real(0.);
        test_res.real(0.);
        test_res.imag(0.);
        for (int x = 0; x < xSize; ++x) {
            bess_sum = 0.;
            for (int i = 0; i <= num_p; ++i) {
                test_ori.imag(v);
                test_res = boost::math::sph_bessel(0, v) * exp(test_ori);
                std::cout << test_res
                          << std::endl;
                v += delta_v;
                bess_sum += bess_tmp * delta_v;
                ori_vec.push_back(bess_tmp);
                int_vec.push_back(bess_sum);
            }
        }
    }
    catch (std::exception ex)
    {
        std::cout << "Thrown exception : " << ex.what() << std::endl;
    }
    return 0;
}

int born_wolf_zstack(std::vector<int> zs,
                     std::vector<std::vector<std::vector<double> > >& M3D,
                     double k, double NA, double n_i, int num_p)
{
    int step = 800 / num_p;
    int z_num = zs.size();
    double bessel_res = 0.0;
    try
    {
        for (int z = 0; z < num_p; z++)
        {
            M3D[k].resize(z_num);
            for (int i = 0; i < num_p; i++)
            {
                M3D[k][i].resize(num_p);
                for (int j = 0; j <= i; j++)
                {
                    bessel_res = born_wolf_point(k, NA, n_i, j*step, i*step, zs[k]);
                    M3D[z][i][j] = bessel_res;
                    M3D[z][j][i] = bessel_res;
                }
            }
            std::cout << " " << z << " / " << z_num << "\tfinished" << std::endl;
        }
    }
    catch (std::exception ex)
    {
        std::cout << "Thrown exception : " << ex.what() << std::endl;
    }
    return 0;
}
