//
// Created by peo on 17-7-13.
//

#include "deconvpsf.h"

int deconv_im(const std::vector<std::vector<double> > raw_im,
              std::vector<std::vector<double> >&dis_im)
{
    int raw_size = raw_im.size();
    int core_size = raw_size * 2;
    dis_im.resize(core_size);
    for (int i = 0; i < core_size; ++i) {
        dis_im[i].resize(core_size);
    }
    fftw_complex *in, *out;
//    fftw_plan_dft_2d(raw_size, raw_size,);
}

int cv_deconv_im(cv::InputArray in_im, cv::OutputArray out_im)
{
    cv::dft(in_im, out_im, cv::DFT_COMPLEX_OUTPUT);
//    std::cout << out_im <<std::endl;
    return 0;
}

fftw_complex* setupHMM(std::vector<std::vector<double> > &vals,
                        int N, int M)
{
    fftw_complex * temp;
    temp = new fftw_complex[N*M];
    for (unsigned i = 0; i < N; ++i) {
        for (unsigned j = 0; j < M; ++j) {
            temp[i*N+j][0] = vals[i][j];
        }
    }
    return temp;
}

void dftplane(int squareSize,
              cv::InputArray in_img,
              cv::OutputArray)
{
    fftw_plan planR, planG, planB;
    fftw_complex *inR, *inG, *inB, *outR, *outG, *outB;

    // allocate input arrays
    inR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    inG = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    inB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);

//    uchar *ptr = in_img.data;
//    for (int i = 0; i < in_img.rows; ++i) {
//        for (int j = 0; j < in_img.cols; ++j) {
//            inB[i*in_img.cols + j][0] = ptr[i*in_img.cols*3 + j*3];
//            inG[i*in_img.cols + j][0] = ptr[i*in_img.cols*3 + j*3+1];
//            inR[i*in_img.cols + j][0] = ptr[i*in_img.cols*3 + j*3+2];
//        }
//    }


    // allocate output arrays
    outR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    outG = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    outB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);

    // create plans
    planR = fftw_plan_dft_2d(squareSize, squareSize, inR, outR, FFTW_FORWARD, FFTW_ESTIMATE);
    planG = fftw_plan_dft_2d(squareSize, squareSize, inG, outG, FFTW_FORWARD, FFTW_ESTIMATE);
    planB = fftw_plan_dft_2d(squareSize, squareSize, inB, outB, FFTW_FORWARD, FFTW_ESTIMATE);

    // perform FORWARD fft
    fftw_execute(planR);
    fftw_execute(planG);
    fftw_execute(planB);
//    std::cout << outR[0][0] << " "
//              << outR[0][1] << std::endl;

    // free memory
    fftw_destroy_plan(planB);
    fftw_destroy_plan(planG);
    fftw_destroy_plan(planR);
    fftw_free(inR); fftw_free(outR);
    fftw_free(inG); fftw_free(outG);
    fftw_free(inB); fftw_free(outB);
}

void iftplane(int squareSize,
             cv::InputArray in_img,
             cv::OutputArray out_img)
{
    fftw_plan planR, planG, planB;
    fftw_complex *inR, *inG, *inB, *outR, *outG, *outB;

    // allocate input arrays
    inR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    inG = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    inB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    // allocate output arrays
    outR = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    outG = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    outB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);

    // create plans
    planR = fftw_plan_dft_2d(squareSize, squareSize, inR, outR, FFTW_FORWARD, FFTW_ESTIMATE);
    planG = fftw_plan_dft_2d(squareSize, squareSize, inG, outG, FFTW_FORWARD, FFTW_ESTIMATE);
    planB = fftw_plan_dft_2d(squareSize, squareSize, inB, outB, FFTW_FORWARD, FFTW_ESTIMATE);

    // perform FORWARD fft
    fftw_execute(planB);
    fftw_execute(planG);
    fftw_execute(planR);

    // free memory
    fftw_destroy_plan(planB);
    fftw_destroy_plan(planG);
    fftw_destroy_plan(planR);
    fftw_free(inR); fftw_free(outR);
    fftw_free(inG); fftw_free(outG);
    fftw_free(inB); fftw_free(outB);
}

void convolutionDFT(cv::InputArray A, cv::InputArray B, cv::OutputArray C)
{
    // A choice whether reallocate the output array
    C.create(abs(A.rows()-B.rows())+1, abs(A.cols() - B.cols())+1, A.type());
    cv::Size dftSize;

    // calculate the size of DFT transform
    dftSize.width = cv::getOptimalDFTSize(A.cols()+B.cols() - 1);
    dftSize.height = cv::getOptimalDFTSize(A.rows()+B.rows() - 1);

    // allocate temporay buffers and initialize them with 0's
    cv::Mat tempA(dftSize, A.type(), cv::Scalar::all(0));
    cv::Mat tempB(dftSize, B.type(), cv::Scalar::all(0));

    // copy A and B to the top-left corner of tempA and tempB, respectively
    cv::Mat roiA(tempA, cv::Rect(0, 0, A.cols(), A.rows()));
    A.copyTo(roiA);
    cv::Mat roiB(tempB, cv::Rect(0, 0, B.cols(), B.rows()));
    B.copyTo(roiB);

    // now transform the padded A & B in-place
    // use "nonzerosRows" hint for faster processing
    cv::dft(tempA, tempA, 0, A.rows());
    cv::dft(tempB, tempB, 0, B.rows());

    // multiply the spectrum
    // the function handles packed spectrum representations well
    cv::mulSpectrums(tempA, tempB, tempA, cv::DFT_COMPLEX_OUTPUT);

    cv::dft(tempA, tempA, cv::DFT_INVERSE + cv::DFT_SCALE, C.rows());

    tempA(cv::Rect(0, 0, C.cols(), C.rows())).copyTo(C);
}

void RichardLucydeconv(cv::InputArray img,
                       cv::InputArray core,
                       cv::OutputArray out_img)
{
    cv::Mat im_correction, im_new_est;
    for (int i = 0; i < 10; ++i) {
        cv::filter2D(img, out_img, img.depth(), core);
        im_correction.create(out_img.size(), out_img.type());
        cv::subtract(img, out_img, im_correction, 0);
        im_new_est.create(out_img.size(), out_img.type());
        cv::add(img, im_correction, im_new_est, 0);
        im_new_est.copyTo(img.getMat());
    }
}