//
// Created by peo on 17-7-13.
//

#include "deconvpsf.h"


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
              cv::OutputArray out_img)
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

void idftplane(int squareSize,
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

int single_fft(cv::InputArray in_img, cv::OutputArray out_img)
{
    fftw_plan plan;
    fftw_complex *inI, *outI;
    int squareSize = in_img.cols() > in_img.rows() ? in_img.cols() : in_img.rows();
    inI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);
    outI = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*squareSize*squareSize);

    for (int i = 0; i < in_img.rows(); i++) {
        for (int j = 0; j < in_img.cols(); j++) {
             outI[i*in_img.rows()+j][0] = in_img.getMat().at<double>(i, j);
             outI[i*in_img.rows()+j][1] = 0.0;
        }
    }

    plan = fftw_plan_dft_2d(squareSize, squareSize, inI, outI, FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan);

    out_img.create(squareSize, squareSize, CV_64F);

    for (int i = 0, k = 0; i < squareSize; i++) {
        for (int j = 0; j < squareSize; j++) {
            out_img.getMat().at<double>(i, j) = outI[k++][0];
            std::cout << outI[k++][0];
        }
        std::cout << std::endl;
    }

    fftw_free(inI);
    fftw_free(outI);

    return 0;
}


void RichardLucy(cv::InputArray img,
                       cv::InputArray core,
                       cv::OutputArray out_img)
{
    cv::Mat im_correction, im_new_est, in_padded, core_padded, adj_padded;
    cv::Mat in_img = img.getMat(), core_img = core.getMat();
    cv::Mat in_complex, core_complex, core_adj, core_conj, norm_img;
    im_new_est = in_img.clone();

    // padded the image for better fft performance
    int width = img.cols() + core.cols();
    int height = img.rows() + core.rows();

    int opt_width = cv::getOptimalDFTSize(width);
    int opt_height = cv::getOptimalDFTSize(height);

    // padded raw image dosen't hurt image performance
    cv::copyMakeBorder(core_img, core_padded, 0, opt_height-core.rows(), 0, opt_width-core.cols(), cv::BORDER_CONSTANT, cv::Scalar::all(0));


    // PSF core padded for computation
    cv::Mat core_planes[] = {cv::Mat_<double>(core_padded), cv::Mat::zeros(core_padded.size(), CV_64F)};
    // PSF core adjoint image plane
    cv::Mat corr_planes[] = {cv::Mat_<double>(in_padded), cv::Mat::zeros(core_padded.size(), CV_64F)};
    cv::Mat in_planes[] = {cv::Mat_<double>(in_padded), cv::Mat::zeros(in_padded.size(), CV_64F)};

    cv::merge(core_planes, 2, core_complex);
    cv::dft(core_complex, core_complex, cv::DFT_COMPLEX_OUTPUT);
    cv::split(core_complex, core_planes);
    core_planes[1] *= -1;
    cv::merge(core_planes, 2, core_conj);
    cv::Rect rect_in = cv::Rect(0, 0, img.cols(), img.rows());

    for (int i = 0; i < 4; i++) {

        // Loop begin
        // Replace the original data
        cv::copyMakeBorder(im_new_est, in_padded, 0, opt_height-img.rows(), 0, opt_width-img.cols(), cv::BORDER_CONSTANT, cv::Scalar::all(0));

        in_planes[0] = in_padded;
        in_planes[1] = cv::Mat::zeros(in_padded.size(), CV_64F);

        // compute the DFT transform
        cv::merge(in_planes, 2, in_complex);
        cv::dft(in_complex, in_complex, cv::DFT_COMPLEX_OUTPUT);
        cv::mulSpectrums(in_complex, core_complex, in_complex, cv::DFT_COMPLEX_OUTPUT);

        // get inverse DFT transform
        cv::idft(in_complex, in_complex);
        cv::split(in_complex, corr_planes);
        cv::magnitude(corr_planes[0], corr_planes[1], corr_planes[0]);
        corr_planes[0](rect_in).copyTo(im_correction);
        cv::divide(in_img, im_correction, im_correction);

        cv::copyMakeBorder(im_correction, in_padded, 0, opt_height-img.rows(), 0, opt_width-img.cols(), cv::BORDER_CONSTANT, cv::Scalar::all(0));
        in_planes[0] = in_padded;
        in_planes[1] = cv::Mat::zeros(in_padded.size(), CV_64F);
        cv::merge(in_planes, 2, in_complex);
        cv::dft(in_complex, in_complex);

        // deconvolution $H^T$
        divSpectrums(in_complex.clone(), core_complex, in_complex, 0, true);
        cv::mulSpectrums(in_complex.clone(), core_complex, in_complex, 0, true);

        cv::idft(in_complex, in_complex);
        cv::split(in_complex, corr_planes);
        cv::magnitude(corr_planes[0], corr_planes[1], corr_planes[0]);
        corr_planes[0](rect_in).copyTo(im_correction);
        cv::multiply(im_new_est, im_correction, im_new_est);
        im_new_est.copyTo(norm_img);
    }
    // compute the convolution in frequency domain:denom
    cv::normalize(norm_img, norm_img, 0, 1, CV_MINMAX);
    norm_img.copyTo(out_img);
}

void fourier_show(cv::InputArray img)
{
    // padded the image for better fft processing
    cv::Mat padded;
    int opt_width = cv::getOptimalDFTSize(img.rows());
    int opt_height = cv::getOptimalDFTSize(img.cols());
    cv::copyMakeBorder(img, padded, 0, opt_width-img.rows(), 0, opt_height-img.cols(), cv::BORDER_CONSTANT, cv::Scalar::all(0));

    // DFT
    cv::Mat planes[] = {cv::Mat_<double>(padded), cv::Mat::zeros(padded.size(), CV_64F)};
    cv::Mat com_img;
    cv::merge(planes, 2, com_img);
    cv::dft(com_img, com_img);

    // get DFT image
    cv::split(com_img, planes);
    cv::magnitude(planes[0], planes[1], planes[0]);

    cv::Mat mag_mat = planes[0];
    mag_mat += cv::Scalar::all(1);
    cv::log(mag_mat, mag_mat);

    // relocation quad image
    mag_mat = mag_mat(cv::Rect(0, 0, mag_mat.cols & -2, mag_mat.rows & -2));
    int cx = mag_mat.cols / 2;
    int cy = mag_mat.rows / 2;

    cv::Mat q0(mag_mat, cv::Rect(0, 0, cx, cy));
    cv::Mat q1(mag_mat, cv::Rect(0, cy, cx, cy));
    cv::Mat q2(mag_mat, cv::Rect(cx, cy, cx, cy));
    cv::Mat q3(mag_mat, cv::Rect(cx, 0, cx, cy));

    cv::Mat tmp;
    q0.copyTo(tmp);
    q2.copyTo(q0);
    tmp.copyTo(q2);

    q1.copyTo(tmp);
    q3.copyTo(q1);
    tmp.copyTo(q3);

    // prepare for exhibition
    cv::Mat mag_img(mag_mat);
    cv::normalize(img.getMat(), img.getMat(), 255, 0);
    cv::normalize(mag_img, mag_img, 255, 0);
    cv::imshow("FOURIER", mag_img);
}

void multiply_fourier(cv::InputArray A, cv::InputArray B, cv::OutputArray C)
{
    cv::Mat paddedA, paddedB, spectrumC;
    cv::Mat matA, matB;
    matA = A.getMat();matB = B.getMat();
    // make more space for preventing image distoration
    int width = A.cols() + B.cols()*2 - 1;
    int height = A.rows() + B.rows()*2 - 1;
    int opt_width = cv::getOptimalDFTSize(width);
    int opt_height = cv::getOptimalDFTSize(height);
    cv::Rect rect(0, 0, matA.cols+matB.cols, matA.rows+matB.rows);

    cv::copyMakeBorder(matA, paddedA, matB.rows, opt_height-matA.rows-matB.rows, matB.cols, opt_width-matA.cols-matB.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::copyMakeBorder(matB, paddedB, opt_height-matB.rows, 0, opt_width-matB.cols, 0, cv::BORDER_CONSTANT, cv::Scalar::all(0));

    // DFT
    cv::Mat planesA[] = {cv::Mat_<double>(paddedA), cv::Mat::zeros(paddedA.size(), CV_64F)};
    cv::Mat planesB[] = {cv::Mat_<double>(paddedB), cv::Mat::zeros(paddedB.size(), CV_64F)};
    cv::Mat planesC[] = {cv::Mat::zeros(paddedB.size(), CV_64F), cv::Mat::zeros(paddedB.size(), CV_64F)};

    cv::Mat com_A, com_B, com_C;
    cv::merge(planesA, 2, com_A);
    cv::merge(planesB, 2, com_B);
    cv::dft(com_A, com_A);
    cv::dft(com_B, com_B);

    cv::mulSpectrums(com_A, com_B, com_C, 0, false);

    cv::dft(com_C, com_C, cv::DFT_INVERSE);
    cv::split(com_C, planesC);
    cv::magnitude(planesC[0], planesC[1], planesC[0]);
    cv::normalize(planesC[0], planesC[0], 0, 1, CV_MINMAX);
    planesC[0](rect).copyTo(C);
//    planesC[0].copyTo(C);
}

void relocation(cv::InputArray A, cv::OutputArray B)
{
    cv::Mat A_mat = A.getMat();
    A_mat = A_mat(cv::Rect(0, 0, A_mat.cols & -2, A_mat.rows & -2));
    int cx = A_mat.cols / 2;
    int cy = A_mat.rows / 2;

    cv::Mat q0(A_mat, cv::Rect(0, 0, cx, cy));
    cv::Mat q1(A_mat, cv::Rect(0, cy, cx, cy));
    cv::Mat q2(A_mat, cv::Rect(cx, cy, cx, cy));
    cv::Mat q3(A_mat, cv::Rect(cx, 0, cx, cy));

    cv::Mat tmp;
    q0.copyTo(tmp);
    q2.copyTo(q0);
    tmp.copyTo(q2);

    q1.copyTo(tmp);
    q3.copyTo(q1);
    tmp.copyTo(q3);
    A_mat.copyTo(B);
}

void divide_fourier(cv::InputArray A, cv::InputArray B, cv::OutputArray C)
{
    cv::Mat paddedA, paddedB, spectrumC, paddedBs;
    cv::Mat matA, matB;
    matA = A.getMat();matB = B.getMat();
    int width = matA.cols;//+matB.cols;
    int height = matA.rows;//+matB.cols;
    int opt_width = width;//cv::getOptimalDFTSize(width);
    int opt_height = height;//cv::getOptimalDFTSize(height);
    cv::Rect rect(0, 0, matA.cols-matB.cols, matA.rows-matB.rows);

    cv::copyMakeBorder(matA, paddedA, 0, opt_height-matA.rows, 0, opt_width-matA.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::copyMakeBorder(matB, paddedB, 0, opt_height-matB.rows, 0, opt_width-matB.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));

    // DFT
    cv::Mat planesA[] = {cv::Mat_<double>(paddedA), cv::Mat::zeros(paddedA.size(), CV_64F)};
    cv::Mat planesB[] = {cv::Mat_<double>(paddedB), cv::Mat::zeros(paddedB.size(), CV_64F)};
    cv::Mat planesC[] = {cv::Mat::zeros(paddedB.size(), CV_64F), cv::Mat::zeros(paddedB.size(), CV_64F)};

    cv::Mat com_A, com_B, com_C, icom_B;
    cv::merge(planesA, 2, com_A);
    cv::merge(planesB, 2, com_B);
    cv::dft(com_A, com_A);
    cv::dft(com_B, com_B);
    cv::split(com_B, planesB);

    divSpectrums(com_A, com_B, com_C, 0, true);

    cv::idft(com_C, com_C);
    cv::split(com_C, planesC);

    cv::magnitude(planesC[0], planesC[1], planesC[0]);
    cv::normalize(planesC[0], planesC[0], 0, 1, CV_MINMAX);
//    planesC[0](rect).copyTo(C);
    planesC[0].copyTo(C);
}

void divide_fourier_test(cv::InputArray A, cv::InputArray B, cv::OutputArray C)
{
    cv::Mat paddedA, paddedB, spectrumC, paddedBs;
    cv::Mat matA, matB;
    matA = A.getMat(), matB = B.getMat();
    int width = matA.cols+matB.cols-1;
    int height = matA.rows+matB.rows-1;
    int opt_width = cv::getOptimalDFTSize(width);
    int opt_height = cv::getOptimalDFTSize(height);
    cv::Rect rect(0, 0, matA.cols-matB.cols, matA.rows-matB.rows);

    cv::copyMakeBorder(matA, paddedA, 0, opt_height-matA.rows, 0, opt_width-matA.cols, cv::BORDER_CONSTANT, cv::Scalar::all(0));
    cv::copyMakeBorder(matB, paddedB, opt_height-matB.rows, 0, opt_width-matB.cols, 0, cv::BORDER_CONSTANT, cv::Scalar::all(0));

    // DFT
    cv::Mat planesA[] = {cv::Mat_<double>(paddedA), cv::Mat::zeros(paddedA.size(), CV_64F)};
    cv::Mat planesB[] = {cv::Mat_<double>(paddedB), cv::Mat::zeros(paddedB.size(), CV_64F)};
    cv::Mat planesC[] = {cv::Mat::zeros(paddedB.size(), CV_64F), cv::Mat::zeros(paddedB.size(), CV_64F)};

    cv::Mat com_A, com_B, com_C, icom_B;
    cv::merge(planesA, 2, com_A);
    cv::merge(planesB, 2, com_B);
    cv::dft(com_A, com_A);
    cv::dft(com_B, com_B);
    cv::split(com_B, planesB);

    divSpectrums(com_A, com_B, com_C, 0, true);

    cv::idft(com_C, com_C);
    cv::split(com_C, planesC);

    cv::magnitude(planesC[0], planesC[1], planesC[0]);
    cv::normalize(planesC[0], planesC[0], 0, 1, CV_MINMAX);
//    planesC[0](rect).copyTo(C);
    planesC[0].copyTo(C);
}

static void divSpectrums(cv::InputArray _srcA, cv::InputArray _srcB, cv::OutputArray _dst, int flags, bool conjB=true)
{
    cv::Mat srcA = _srcA.getMat(), srcB = _srcB.getMat();
    int depth = srcA.depth(), cn = srcA.channels(), type = srcA.type();
    int rows = srcA.rows, cols = srcA.cols;

    CV_Assert( type == srcB.type() && srcA.size() == srcB.size() );
    CV_Assert( type == CV_64FC1 || type == CV_64FC2);

    _dst.create(srcA.rows, srcA.cols, type);
    cv::Mat dst = _dst.getMat();

    CV_Assert(dst.data != srcA.data); // non-inplace check
    CV_Assert(dst.data != srcB.data); // non-inplace check

    bool is_1d = (flags & cv::DFT_ROWS) || (rows == 1 || (cols == 1 &&
            srcA.isContinuous() && srcB.isContinuous() && dst.isContinuous()));

    if (is_1d && !(flags & cv::DFT_ROWS)) {
        cols = cols + rows - 1, rows = 1;
    }
    int ncols = cols*cn;
    int j0 = cn == 1;
    int j1 = ncols - (cols % 2 == 0 && cn == 1);

    if ( depth == CV_64F ) {
        const double *dataA = srcA.ptr<double>();
        const double *dataB = srcB.ptr<double>();
        auto *dataC = dst.ptr<double>();
        auto eps = DBL_EPSILON;

        size_t stepA = srcA.step / sizeof(dataA[0]);
        size_t stepB = srcB.step / sizeof(dataB[0]);
        size_t stepC = dst.step / sizeof(dataC[0]);

        if (!is_1d && cn == 1) {
            // two channels real and image
            for (int k = 0; k < (cols % 2 ? 1 : 2); k++) {
                if (k == 1) {
                    dataA += cols - 1, dataB += cols - 1, dataC += cols - 1;
                }
                dataC[0] = dataA[0] / (dataB[0] + eps);
                if (rows % 2 == 0) {
                    dataC[(rows - 1) * stepC] = dataA[(rows - 1) * stepA] / (dataB[(rows - 1) * stepB] + eps);
                }
                if (!conjB) {
                    for (int j = 1; j <= rows - 2; j += 2) {
                        double denom =
                                dataB[j * stepB] * dataB[j * stepB] + dataB[(j + 1) * stepB] * dataB[(j + 1) * stepB] +
                                eps;
                        double re =
                                dataA[j * stepA] * dataB[j * stepB] + dataA[(j + 1) * stepA] * dataB[(j + 1) * stepB];
                        double im =
                                dataA[(j + 1) * stepA] * dataB[j * stepB] - dataA[j * stepA] * dataB[(j + 1) * stepB];

                        dataC[j * stepC] = (re / denom);
                        dataC[(j + 1) * stepC] = (im / denom);
                    }
                } else {
                    for (int j = 1; j <= rows - 2; j += 2) {
                        double denom =
                                dataB[j * stepB] * dataB[j * stepB] + dataB[(j + 1) * stepB] * dataB[(j + 1) * stepB] +
                                eps;
                        double re =
                                dataA[j * stepA] * dataB[j * stepB] - dataA[(j + 1) * stepA] * dataB[(j + 1) * stepB];
                        double im =
                                dataA[(j + 1) * stepA] * dataB[j * stepB] + dataA[j * stepA] * dataB[(j + 1) * stepB];

                        dataC[j * stepC] = re / denom;
                        dataC[(j + 1) * stepC] = im / denom;
                    }
                }
                if (k == 1) {
                    dataA -= cols - 1, dataB -= cols - 1, dataC -= cols - 1;
                }
            }
        }

        for (; rows--; dataA += stepA, dataB += stepB, dataC += stepC) {
            if (is_1d && cn == 1) {
                dataC[0] = dataA[0] / (dataB[0] + eps);
                if (cols % 2 == 0) {
                    dataC[j1] = dataA[j1] / (dataB[j1] + eps);
                }
            }
            if (!conjB) {
                for (int j = 0; j < j1; j += 2) {
                    double denom = dataB[j] * dataB[j] + dataB[j + 1] * dataB[j + 1] + eps;
                    double re = dataA[j] * dataB[j] - dataA[j + 1] * dataB[j + 1];
                    double im = dataA[j + 1] * dataB[j] + dataA[j] * dataB[j + 1];
                    dataC[j] = re / denom;
                    dataC[j + 1] = im / denom;
                }
            } else {
                for (int j = j0; j < j1; j += 2) {
                    double denom = dataB[j] * dataB[j] + dataB[j + 1] * dataB[j + 1] + eps;
                    double re = dataA[j] * dataB[j] + dataA[j + 1] * dataB[j + 1];
                    double im = dataA[j + 1] * dataB[j] - dataA[j] * dataB[j + 1];
                    dataC[j] = re / denom;
                    dataC[j + 1] = im / denom;
                }
            }
        }
    }
}