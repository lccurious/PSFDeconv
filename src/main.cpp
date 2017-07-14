//
// Created by peo on 17-7-5.
//

#include <boost/program_options.hpp>
#include <boost/program_options/errors.hpp>
#include <string>
#include <iostream>

#include "cvutils.h"
#include "psf_gen2.h"
#include "comparemat.h"
#include "vtkutils.h"
#include "deconvpsf.h"
#include "config.h"

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

#if defined(PSF_DEMON) || defined(INTEGRAL_TEST) || defined(PSF_GEN_COMPARE)// show psf demon
    std::vector<std::vector<double> > psf_matrix;
#endif // PSF_DEMON

#ifdef INTEGRAL_TEST
    std::vector<double> ori_vec, int_vec;
    integral_test(0, ori_vec, int_vec, 1000, 56);
    vtk_2Dplot_com(ori_vec, int_vec);
#endif // integral TEST PART

#ifdef PSF_DEMON
    boost::progress_display *show_progress = NULL;
    show_progress = new boost::progress_display(stack_depth);
    born_wolf_full((stack_depth-32)*170, psf_matrix, M_2PI/em_wavelen, NA, refr_index, psf_size);

    for(int i = 0; i < stack_depth; i++) {
        ++(*show_progress);
    }
#endif // show psf demon end

#ifdef PSF_GEN_COMPARE
    std::vector<double> X;
    std::vector<std::vector<double> > M2D;
    born_wolf_full((stack_depth-32)*270, psf_matrix, M_2PI/em_wavelen, NA, refr_index, psf_size);
    mat2vector(compare_filename.c_str(), M2D, stack_depth);
    for (int i = 0; i < psf_matrix[0].size(); i++) {
        X.push_back((double)i);
    }

#if defined(PSF_GEN_COMPARE) && defined(_PRINT_M2D)
    for (int i = 0; i < M2D.size(); i++) {
        for (int j = 0; j < M2D[0].size(); j++) {
            std::cout << M2D[i][j] << "\t";
        }
        std::cout << std::endl;
    }
#endif

    vtk_2Dplot_com(psf_matrix[127], M2D[127]);
//    vtk_2Dplot(X, M2D[127]);
//    simple_show(psf_matrix);
#endif // PSF_GEN_COMPARE

#ifdef DFT_TEST
    cv::Mat in_img, out_img;
    in_img = cv::imread("images/B0008393-800.jpg", 0);
    int squareSize = in_img.cols < in_img.rows ? in_img.rows:in_img.cols;
    cv::Mat padded;

    // pad the raw image can make the transformation more effective
    int m = cv::getOptimalDFTSize(in_img.rows);
    int n = cv::getOptimalDFTSize(in_img.cols);
    cv::copyMakeBorder(in_img, padded, 0, m-in_img.rows, 0, n-in_img.cols,
                       cv::BORDER_CONSTANT, cv::Scalar::all(0));

    std::cout << "m = " << m << " rows = " << in_img.rows
              << "n = " << n << " cols = " << in_img.cols
              << std::endl;
    // planes with REAL part and the IMAGE part
    cv::Mat planes[] = {cv::Mat_<double>(padded), cv::Mat::zeros(padded.size(), CV_64F)};
    cv::Mat complexI;
    std::vector<cv::Mat> mv;
    cv::merge(planes, 2, complexI);
    cv::dft(complexI, complexI);
    cv::idft(complexI, padded);
    cv::split(complexI, mv);
    cv::Mat manplane;
    cv::magnitude(mv[0], mv[1], manplane);
    cv::imshow("DFT", manplane);
    cv::waitKey(-1);
    std::cout << "assigned successful!" << std::endl;
#endif // DFT TEST

#ifdef LOAD_TIFF_TEST
    cv::Mat image;
    image = cv::imread("/media/peo/Docunment/DeskTop/deconv/SmallFOV/Actin_CoverslipView.tif");
    if (! image.data)
    {
        std::cout << "Could not open or find the image" << std::endl;
        return -1;
    }
    cv::namedWindow("TIFF", CV_WINDOW_AUTOSIZE);
    cv::imshow("TIFF", image);
    cv::waitKey(-1);
#endif

#ifdef DECONVOLUTION_TEST

#ifndef STACK_SIZE_RATIO
#define STACK_SIZE_RATIO 200
#endif

    cv::Mat in_image, out_image, psf_core;
    in_image = cv::imread("images/sdaefawefawefceacaec.png", 0);
    std::vector<std::vector<double> > psf_matrix;
    born_wolf_full((stack_depth-32)*STACK_SIZE_RATIO, psf_matrix, M_2PI/em_wavelen,NA,refr_index,psf_size);
    psf_core.create(psf_matrix.size(), psf_matrix.size(), CV_64FC2);
    std::cout << "PSF Size: " << psf_matrix.size()
              << "\nCORE SIZE: " << psf_core.rows << "\t" << psf_core.cols
              << std::endl;
    vec2mat(psf_matrix, psf_core);
    in_image.convertTo(in_image, CV_64FC2);
    std::cout << "transformed channel: " << in_image.channels() << std::endl;
    cv::imshow("RAW_IMG", in_image);
    convolutionDFT(in_image, psf_core, out_image);
    cv::imshow("PSF_CORE", out_image);
    cv::waitKey(-1);
#endif

#ifdef TIFF_READ_TEST
    cv::Mat TIFF_img;
    getTIFF("/media/peo/Docunment/DeskTop/PSF/deconv/SmallFOV/Actin_CoverslipView.tif", TIFF_img);
#endif

    return 0;
}
