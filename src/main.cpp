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
    int step = AiryRadius / psf_size;
    born_wolf_full((stack_depth-32)*step, psf_matrix, M_2PI/em_wavelen, NA, refr_index, psf_size);
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

    vtk_2Dplot_com(psf_matrix[127], M2D[127], "Local", "ImageJ");
//    vtk_2Dplot(X, M2D[127]);
    simple_show_vec(psf_matrix);
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

    cv::Mat in_image, out_image, psf_core, psf_core_inv,complex_img, psf_show;
    std::vector<std::vector<double> > psf_matrix;

    // prepare the psf core for convolution
    born_wolf_full((stack_depth-32), psf_matrix, M_2PI/em_wavelen,NA,refr_index,psf_size);
    psf_core.create(psf_matrix.size(), psf_matrix.size(), CV_64F);
    vec2mat(psf_matrix, psf_core);
    cv::normalize(psf_core, psf_core, 0, 1, CV_MINMAX);
    cv::resize(psf_core, psf_show, cv::Size(512, 512));

    // Read raw TIFF format
    int total_seq = TIFFframenumber(test_image_name.c_str());
    getTIFF(test_image_name.c_str(), in_image, 0);
    int width = in_image.cols, height = in_image.rows;

    // for fitting the screen size
    while (width > 1000 || height > 1000) {
        width >>= 1;
        height >>= 1;
    }
    int k = 0;
    char image_name[50];
    boost::progress_display *show_progress = NULL;
    show_progress = new boost::progress_display(total_seq);
    for (int i = 0; i < total_seq; ++i) {
        getTIFF(test_image_name.c_str(), in_image, i);

        // applying the PSF convolution operation
        // TODO(peo):Set up denoise mechanism
//         divide_fourier(in_image, psf_core, out_image);
//        RichardLucydeconv(in_image.clone(), psf_core, out_image);
        RichardLucy(in_image.clone(), psf_core.clone(), out_image);

        // operation for exhibition
        cv::normalize(in_image, in_image, 0, 1, CV_MINMAX);
        cv::resize(in_image, in_image, cv::Size(width, height));
        cv::resize(out_image, out_image, cv::Size(width, height));
        cv::imshow("RAW", in_image);
        cv::imshow("PSF", psf_show);
        cv::imshow("DE", out_image);
        k = cv::waitKey(60);
        if (k == 27) {
            std::cout << std::endl;
            break;
        }
        if (k == ' ') {
            std::sprintf(image_name, "logs/Ricard-Lucy_%d.png", i);
            out_image *= 255;
            out_image.convertTo(out_image, CV_8UC1);
            cv::imwrite(image_name, out_image);
        }
        ++(*show_progress);
    }
#endif

#ifdef TIFF_READ_TEST
    getTIFFinfo(test_image_name.c_str());
    cv::Mat TIFF_img;
    int k = 0;
    int TIFF_frame_num = TIFFframenumber(test_image_name.c_str());
    for (int i = 0; i < TIFF_frame_num; ++i) {
        getTIFF(test_image_name.c_str(), TIFF_img, i);
        cv::normalize(TIFF_img, TIFF_img, 255, 0);
        cv::putText(TIFF_img, std::to_string(i).c_str(), cv::Point(25, 55),
                    cv::FONT_HERSHEY_DUPLEX, 2.0, CV_RGB(255, 255, 255), 2);
        cv::resize(TIFF_img, TIFF_img, cv::Size(TIFF_img.cols/2, TIFF_img.rows/2),
                   0.5, 0.5);
        cv::imshow("TIFF Frame", TIFF_img);
        k = cv::waitKey(25);
        if (k == 27) {
            break;
        }
    }
#endif

#ifdef TRACK_BAR_NEURON
    cv_slideWin(test_image_name.c_str());
#endif

#ifdef PSF_STACK_VIEW
    cv::Mat psf_core, stack_view;
    std::vector<std::vector<double> > psf_matrix;
    double step = AiryRadius / psf_size;
    double bessel_res = 0.0;
    stack_view.create(psf_size*2, 64, CV_64F);
    for (int i = 0; i < 64; i++) {
        for (int j = 0; j < psf_size; j++) {
            bessel_res = born_wolf_point(M_2PI/em_wavelen, NA, refr_index, j*step, 0, (i-32)*step);
            stack_view.at<double>(psf_size-j, i) = bessel_res;
            stack_view.at<double>(psf_size+j, i) = bessel_res;
        }
    }
    cv::normalize(stack_view, stack_view, 0, 1, CV_MINMAX);
    stack_view *= 255;
    stack_view.convertTo(stack_view, CV_8UC1);
    cv::resize(stack_view, stack_view, stack_view.size()*2);
    cv::imwrite("logs/psf_stack.png", stack_view);
    cv::imshow("PSF_STACK", stack_view);
    cv::waitKey(-1);
#endif

#ifdef RLAL_TEST
    cv::Mat in_image, out_image, psf_core;
    cv::Mat in_show, out_show, frame, im_gray;
    std::vector<std::vector<double> > psf_matrix;

    born_wolf_full((stack_depth-32), psf_matrix, M_2PI/em_wavelen,NA,refr_index,psf_size);
    psf_core.create((int)psf_matrix.size(), (int)psf_matrix.size(), CV_64F);
    vec2mat(psf_matrix, psf_core);
    double core_sum = cv::sum(psf_core)[0];
    psf_core /= core_sum;

    cv::VideoCapture cap = cv::VideoCapture("/media/peo/Docunment/视频/vlc-record-2017-05-27-19h49m41s-DJI_0012.MP4-.mp4");

    int k, count = (int)cap.get(cv::CAP_PROP_FRAME_COUNT);
    cap >> in_image;
    count --;
    int width = in_image.cols, height = in_image.rows;

    // for fitting the screen size
    while (width > 1000 || height > 1000) {
        width >>= 2;
        height >>= 2;
    }

    // main loop begin
    for(int i = 0; i < count; i++) {
        cap >> frame;
        cv::cvtColor(frame, im_gray, cv::COLOR_BGR2GRAY);
        cv::resize(im_gray, im_gray, cv::Size(width, height));
        multiply_fourier(im_gray, psf_core, in_image);

        // TODO(peo):substitute the function.
        // RichardLucy(in_image, psf_core, out_image);

        divide_fourier(in_image, psf_core, out_image);
        cv::imshow("CONV", in_image);
        cv::imshow("DE", out_image);
        k = cv::waitKey(1);
        if (k == 27) {
            break;
        }
    } // main loop end
    cap.release();
#endif

    return 0;
}
