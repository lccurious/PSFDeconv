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
    int num_p = 10;
    double delta_v = 1.0 / num_p;
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
    }
    catch ( std::exception ex)
    {
        std::cout <<  "Thrown exception " << ex.what() << std::endl;
    }

    return std::abs(bess_sum);
}

/*
double born_wol_acc(std::vector<int>z_stack ,
                    std::vector<std::vector<std::vector<double> > >&M3D,
                    double k, double NA, double n_i, int num_p)
{
    int step = 800 / num_p;
    double bessel_res = 0.0, const_ratio;
    // TODO(peo): Here using magic number for bessel integral accuracy indicator
    double delta_v = 1.0 / 10;
    double x, y;
    std::vector<std::vector<std::vector<double> > >bess_tmp;
    try {
        for (int i = 0; i < num_p; i++) {
            const_ratio.resize(num_p);
            x = i*step;
            for (int j = 0; j < num_p; i++) {
                y = j*step;
                const_ratio = k*NA / n_i * std::sqrt(x*x + y*y);
                for (int v = 0; v < 10; v++) {
                    bess_tmp[i][j].push_back(boost::math::sph_bessel(0, ));
                }
            }
        }
    }
    catch (std::exception ex) {
        std::cout << "Thrown exception: " << ex.what() << std::endl;
    }

}
 */

int born_wolf(int z, std::vector<std::vector<double> >& M2D,
              double k, double NA, double n_i, int num_p)
{
    int step = 800 / num_p;
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
