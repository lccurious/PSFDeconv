//
// Created by peo on 17-6-24.
//
// This file will try to use some multiprocess method to
// Make program faster

#include <boost/thread.hpp>
#include "psf_gen2.h"
#include <boost/progress.hpp>
#include <iostream>

#define M_2PI (6.283185307179586476925286766559)
double ex_wavelen = 488;
double NA = 1.2;
double refr_index = 1.333;
int psf_count = 0;
int psf_pid = 0;
int psf_status;

boost::progress_display *show_progress = NULL;


int psf_step(int z, std::vector<std::vector<double> >& M2D,
          double k, double NA, double n_i, int num_p)
{
    born_wolf(z, M2D, k, NA, n_i, num_p);
    ++(*show_progress);
    // std::cout << "PSF: " << psf_count++ << std::endl;
    psf_status++;
    return 0;
}

int main()
{
    int stack_depth = 512;
    int num_p = 256;
    boost::thread psf_thrd[stack_depth];
    std::vector<std::vector<std::vector<double> > >psf_matrix(stack_depth);
    // Initial the 3D matrix
    for (int i = 0; i < stack_depth; ++i) {
        psf_matrix[i].resize(num_p);
        for (int j = 0; j < num_p; ++j) {
            psf_matrix[i][j].resize(num_p);
        }
    }
    show_progress = new boost::progress_display(stack_depth);

    // Build the thread pool
    clock_t start = clock();
    for (int i = 0; i < stack_depth; ++i) {
        psf_thrd[i] = boost::thread(&psf_step, i, psf_matrix[i],
                                    M_2PI/ex_wavelen, NA, refr_index, num_p);
        // psf_thrd[i].join();
        // psf_step(i, psf_matrix[i], M_2PI/ex_wavelen, NA, refr_index, num_p);
        // Time spent: 116 seconds
        // stack_depth = 200;
    }
    clock_t midum = clock();

    for (int i = 0; i < stack_depth; ++i) {
        psf_thrd[i].join();
    }
    // bool all_done = false;
    // int z = 10;

    /*
    while (!all_done) {
        // std::cout << "PID " << psf_pid << " Joinable "
        //           << psf_thrd[psf_pid].joinable() <<std::endl;
        if (!psf_status[psf_pid]) {
            // std::cout << "PID Joinable" << std::endl;
            psf_thrd[psf_pid] = boost::thread(&psf_step, z, psf_matrix[z],
                                              M_2PI/ex_wavelen, NA, refr_index, num_p, psf_pid);
            psf_thrd[psf_pid].join();
            z += 1;
        }
        psf_pid = (psf_pid + 1) % 10;
        if (z == stack_depth - 1) {
            all_done = true;
        }
    }
     */
    while(psf_status < stack_depth) {}
    clock_t finish = clock();
    std::cout << "Time spent: " << (finish-start) / CLOCKS_PER_SEC
              << "\tBuilt Time: " << (midum-start) / CLOCKS_PER_SEC
              << std::endl;

    return 0;
}