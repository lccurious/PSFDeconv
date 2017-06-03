#include <Eigen/Eigen>
#include <opencv2/opencv.hpp>
#include <vector>
#include <algorithm>
#include <iterator>
#include <opencv2/core.hpp>

using namespace std;
using namespace cv;

void Recomb(Mat &src, Mat &dst)
{
	int cx = src.cols >> 1;
	int cy = src.rows >> 1;
	Mat tmp;
	tmp.create(src.size(), src.type());
	src(Rect(0, 0, cx, cy)).copyTo(tmp(Rect(cx, cy, cx, cy)));
	src(Rect(cx, cy, cx, cy)).copyTo(tmp(Rect(0, 0, cx, cy)));
	src(Rect(cx, 0, cx, cy)).copyTo(tmp(Rect(0, cy, cx, cy)));
	src(Rect(0, cy, cx, cy)).copyTo(tmp(Rect(cx, 0, cx, cy)));
	dst = tmp;
}


void convDFT(Mat& A, Mat& B, Mat& C)
{
	// reallocate the output array if needed
    C.create(abs(A.rows - B.rows) + 1, abs(A.cols - B.cols) + 1, A.type());
    Size dftSize;
    // compute the size of DFT transform
    dftSize.width = getOptimalDFTSize(A.cols + B.cols - 1);
    dftSize.height = getOptimalDFTSize(A.rows + B.rows - 1);
    // allocate temporary buffers and initialize them with 0's
    Mat tempA(dftSize, A.type(), Scalar::all(0));
    Mat tempB(dftSize, B.type(), Scalar::all(0));
    // copy A and B to the top-left corners of tempA and tempB, respectively
    Mat roiA(tempA, Rect(0, 0, A.cols, A.rows));
    A.copyTo(roiA);
    Mat roiB(tempB, Rect(0, 0, B.cols, B.rows));
    B.copyTo(roiB);
    // now transform the padded A & B in-place;
    // use "nonzeroRows" hint for faster processing
    dft(tempA, tempA, 0, A.rows);

    dft(tempB, tempB, 0, A.rows);

    // multiply the spectrums;
    // the function handles packed spectrum representations well
    mulSpectrums(tempA, tempB, tempA, 0);
    // transform the product back from the frequency domain.
    // Even though all the result rows will be non-zero,
    // you need only the first C.rows of them, and thus you
    // pass nonzeroRows == C.rows

    dft(tempA, tempA, DFT_INVERSE + DFT_SCALE);
    // now copy the result back to C.
    C = tempA(Rect((dftSize.width - A.cols) / 2, (dftSize.height - A.rows) / 2, A.cols, A.rows)).clone();
    // all the temporary buffers will be deallocated automatically
}

//----------------------------------------------------------
// Compute Re and Im planes of FFT from Image 
//----------------------------------------------------------
void ForwardFFT(Mat &Src, Mat *FImg)
{
    int M = getOptimalDFTSize(Src.rows);
    int N = getOptimalDFTSize(Src.cols);
    Mat padded;
    copyMakeBorder(Src, padded, 0, M - Src.rows, 0, N - Src.cols, BORDER_CONSTANT, Scalar::all(0));
    Mat planes[] = { Mat_<double>(padded), Mat::zeros(padded.size(), CV_64FC1) };
    Mat complexImg;
    merge(planes, 2, complexImg);
    dft(complexImg, complexImg);
    split(complexImg, planes);
    // crop result
    planes[0] = planes[0](Rect(0, 0, Src.cols, Src.rows));
    planes[1] = planes[1](Rect(0, 0, Src.cols, Src.rows));
    FImg[0] = planes[0].clone();
    FImg[1] = planes[1].clone();
}

//----------------------------------------------------------
// Compute image from Re and Im parts of FFT
//----------------------------------------------------------
void InverseFFT(Mat *FImg, Mat &Dst)
{
    Mat complexImg;
    merge(FImg, 2, complexImg);
    dft(complexImg, complexImg, DFT_INVERSE + DFT_SCALE);
    split(complexImg, FImg);
    Dst = FImg[0];
}

void wienerFilter(Mat &src, Mat &dst, Mat &_h, double k)
{
    //---------------------------------------------------
    // Small epsilon to avoid division by zero
    //---------------------------------------------------
    const double eps = 1E-8;
    //---------------------------------------------------
    int ImgW = src.size().width;
    int ImgH = src.size().height;
    //--------------------------------------------------
    Mat Yf[2];
    ForwardFFT(src, Yf);
    //--------------------------------------------------
    Mat h = Mat::zeros(ImgH, ImgW, CV_64FC1);

    int padx = h.cols - _h.cols;
    int pady = h.rows - _h.rows;

    copyMakeBorder(_h, h, pady / 2, pady - pady / 2, padx / 2, padx - padx / 2, BORDER_CONSTANT, Scalar::all(0));

    Mat Hf[2];
    ForwardFFT(h, Hf);


    //--------------------------------------------------
    Mat Fu[2];
    Fu[0] = Mat::zeros(ImgH, ImgW, CV_64FC1);
    Fu[1] = Mat::zeros(ImgH, ImgW, CV_64FC1);

    complex<double> a;
    complex<double> b;
    complex<double> c;

    double Hf_Re;
    double Hf_Im;
    double Phf;
    double hfz;
    double hz;
    double A;

    for (int i = 0; i < h.rows; i++)
    {
        for (int j = 0; j < h.cols; j++)
        {
            Hf_Re = Hf[0].at<double>(i, j);
            Hf_Im = Hf[1].at<double>(i, j);
            Phf = Hf_Re*Hf_Re + Hf_Im*Hf_Im;
            hfz = (Phf < eps)*eps;
            hz = (h.at<double>(i, j) > 0);
            A = Phf / (Phf + hz + k);
            a = complex<double>(Yf[0].at<double>(i, j), Yf[1].at<double>(i, j));
            b = complex<double>(Hf_Re + hfz, Hf_Im + hfz);
            c = a / b; // Deconvolution :) other work to avoid division by zero
            Fu[0].at<double>(i, j) = (c.real()*A);
            Fu[1].at<double>(i, j) = (c.imag()*A);
        }
    }
    InverseFFT(Fu, dst);
    Recomb(dst, dst);
}

