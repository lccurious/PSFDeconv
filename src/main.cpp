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
    std::vector<std::vector<double> >psf_matrix(psf_size);
    boost::progress_display *show_progress = NULL;
    show_progress = new boost::progress_display(stack_depth);

    for(int i = 0; i < stack_depth; i++) {
        born_wolf(i, psf_matrix, M_2PI/ex_wavelen, NA, refr_index, psf_size);
        ++(*show_progress);
    }
//
    show2DVec(psf_matrix);
    // show psf demon end

//    std::vector<std::vector<std::vector<double> > > M3D;
    std::vector<double> X, Y;
//    mat3vector(compare_filename.c_str(), M3D);
    for (int i = 0; i < psf_matrix[0].size(); i++) {
        X.push_back((double)i);
    }
//    for (int j = 0; j < 256; ++j) {
    vtk_2Dplot(X, psf_matrix[0]);
//    }

    return 0;
}
