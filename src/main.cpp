//
// Created by peo on 17-7-5.
//

#include <boost/program_options.hpp>
#include <boost/program_options/errors.hpp>
#include <string>
#include <iostream>

#include "cvTools.h"
#include "psf_gen2.h"
#include "comparemat.h"
#include "vtkTools.h"

#define M_2PI (6.283185307179586476925286766559)

namespace opt = boost::program_options;

int main(int argc, const char *argv[])
{
    opt::options_description desc("All options");
    desc.add_options()
            ("psf_size,p", opt::value<int>()->default_value(256),
            "the size of PSF kernel")
            ("deconvolution,d", opt::value<std::string>(),
            "the full path of file")
            ("test", opt::value<int>()->default_value(65),
             "test the psf generate cost")
            ("ex_wavelen", opt::value<double>()->default_value(488.0))
            ("em_wavelen", opt::value<double>()->default_value(520.0))
            ("pinhole_radius", opt::value<double>()->default_value(0.55))
            ("refr_index", opt::value<double>()->default_value(1.333))
            ("NA", opt::value<double>()->default_value(1.2))
            ("stack_depth", opt::value<int>()->default_value(65))
            ("compare", opt::value<std::string>(),
            "file contain the ground truth data of BW model")
            ("help", "produce help message")
    ;
    opt::variables_map vm;
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }

    // read parameter from file
    try {
        opt::store(
                opt::parse_config_file<char>("psf.cfg", desc),
                vm
        );
    } catch (const opt::reading_file& e) {
        std::cout << "Failed to open file 'psf.cfg': "
                  << e.what() << std::endl;
    }
    // 解析命令行选项并把值存储到'vm'
    opt::store(opt::parse_command_line(argc, argv, desc), vm);
    opt::notify(vm);

    // the width and height of psf
    int psf_size;
    int stack_depth = vm["stack_depth"].as<int>();
    double NA = vm["NA"].as<double>();
    double refr_index = vm["refr_index"].as<double>();
    double ex_wavelen = vm["ex_wavelen"].as<double>();
    double em_wavelen = vm["em_wavelen"].as<double>();
    double pinhole_radius = vm["pinhole_radius"].as<double>();
    std::string compare_filename = vm["compare"].as<std::string>();

    if (vm.count("psf_size")) {
        psf_size = vm["psf_size"].as<int>();
        std::cout << "start generate "
                  << psf_size << "x" << psf_size
                  << " kernel" << std::endl;
    }

    // show psf demon
    std::vector<std::vector<double> > psf_matrix;
    std::vector<double> ori_vec, int_vec;

    // TODO(peo): integral Test Part
//    integral_test(0, ori_vec, int_vec, 1000, 56);
//    vtk_2Dplot_com(ori_vec, int_vec);
    // integral TEST PART

//    boost::progress_display *show_progress = NULL;
//    show_progress = new boost::progress_display(stack_depth);

//    for(int i = 0; i < stack_depth; i++) {
    born_wolf_full(stack_depth*150, psf_matrix, M_2PI/ex_wavelen, NA, refr_index, psf_size);
//        ++(*show_progress);
//    }
    // show psf demon end

//    std::vector<double> X;
    std::vector<std::vector<double> > M2D;
    mat2vector(compare_filename.c_str(), M2D, stack_depth);
//    for (int i = 0; i < psf_matrix[0].size(); i++) {
//        X.push_back((double)i);
//    }

//    for (int i = 0; i < M2D.size(); i++) {
//        for (int j = 0; j < M2D[0].size(); j++) {
//            std::cout << M2D[i][j] << "\t";
//        }
//        std::cout << std::endl;
//    }
    vtk_2Dplot_com(psf_matrix[127], M2D[127]);
//    vtk_2Dplot(X, M2D[127]);
    simple_show(psf_matrix);
//    vtk_compare_linear(compare_filename.c_str(), psf_matrix[0], M2D[0]);

    return 0;
}
