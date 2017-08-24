//
// Created by peo on 17-7-5.
//

#include <boost/program_options.hpp>
#include <boost/program_options/errors.hpp>
#include <string>
#include <iostream>

#include "cvutils.h"
#include "psfgen.h"
#include "comparemat.h"
#include "vtkutils.h"
#include "deconvpsf.h"
#include "config.h"

#define M_2PI (6.283185307179586476925286766559)

namespace opt = boost::program_options;

int AiryRadius;

int main(int argc, const char *argv[])
{
    opt::options_description desc("All options");
    desc.add_options()
            ("psf_size,p", opt::value<int>()->default_value(256),
            "the size of PSF kernel")
            ("deconvolution,d", opt::value<std::string>(),
            "the full path of file")
            ("test", opt::value<int>()->default_value(32),
             "test the psf generate cost")
            ("ex_wavelen", opt::value<double>()->default_value(488.0))
            ("em_wavelen", opt::value<double>()->default_value(520.0))
            ("pinhole_radius", opt::value<double>()->default_value(0.55))
            ("refr_index", opt::value<double>()->default_value(1.333))
            ("NA", opt::value<double>()->default_value(1.2))
            ("stack_depth", opt::value<int>()->default_value(65))
            ("compare", opt::value<std::string>(),
            "file contain the ground truth data of BW model")
            ("test_image", opt::value<std::string>(),
            "read an image for test algorithm")
            ("logs", opt::value<std::string>()->default_value("logs"),
            "directory for storing generate images")
            ("AiryRadius", opt::value<int>()->default_value(1000))
            ("help", "produce help message")
    ;
    opt::variables_map vm;

    // read parameter from file
    try {
        opt::store(
                opt::parse_config_file<char>("/media/peo/Curious/Visual/deconv/deconv_open/psf.cfg", desc),
                vm
        );
    } catch (const opt::reading_file& e) {
        std::cout << "Failed to open file 'psf.cfg': "
                  << e.what() << std::endl;
    }

    // 解析命令行选项并把值存储到'vm'
    opt::store(opt::parse_command_line(argc, argv, desc), vm, true);

    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 1;
    }
    opt::notify(vm);

    // the width and height of psf
    int psf_size = vm["psf_size"].as<int>();;
    int stack_depth = vm["stack_depth"].as<int>();
    double NA = vm["NA"].as<double>();
    double refr_index = vm["refr_index"].as<double>();
    double ex_wavelen = vm["ex_wavelen"].as<double>();
    double em_wavelen = vm["em_wavelen"].as<double>();
    double pinhole_radius = vm["pinhole_radius"].as<double>();
    AiryRadius = vm["AiryRadius"].as<int>();
    std::string compare_filename = vm["compare"].as<std::string>();
    std::string test_image_name = vm["test_image"].as<std::string>();

    if (vm.count("psf_size")) {
        psf_size = vm["psf_size"].as<int>();
        std::cout << "start generate "
                  << psf_size*2 << "x" << psf_size*2
                  << " kernel" << std::endl;
    }


#ifdef RESTORE_TEST
    cv::Mat in_image, out_image, PSF_tmp;
    std::vector<std::vector<double> > psf_matrix;
    std::vector<std::vector<double> > M2D;
    int k;
    int num_stackin = TIFFframenumber(test_image_name.c_str());
    std::cout << "Start processing " << num_stackin << " TIFF raw image" << std::endl;
    std::cout << "Start generating PSF core..." << std::endl;
    mat2vector(compare_filename.c_str(), M2D, stack_depth);
    vec2mat(M2D, PSF_tmp);
    char RLName[60];
    double *buf = (double *) _TIFFmalloc(in_image.cols * in_image.rows * sizeof(double));

    boost::progress_display *show_progress = NULL;
    show_progress = new boost::progress_display(num_stackin);
    for (int i = 0; i < num_stackin; ++i) {
        getTIFF(test_image_name.c_str(), in_image, i);

        RichardLucy(in_image, PSF_tmp, out_image, 5);
        fourier_show(in_image, "RAWFOURIER");
        norm_show(in_image, "RAW");
        norm_show(out_image, "RL");
        k = cv::waitKey(1);
        if (k == 27) {
            break;
        }
        ++(*show_progress);
    }

#endif

#ifdef TIFF_WRITE_TEST
    cv::Mat image;
    image = cv::imread("images/born_wolf.png", 0);
    cv::imshow("Born", image);
    SingleMatToTiff(image, "born.tiff");
    cv::waitKey(-1);
#endif

    return 0;
}
