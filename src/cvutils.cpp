//
// Created by peo on 17-7-10.
//

#include "cvutils.h"
char *filename_TIFF;
int alpha_slider;


int cv_draw_pie(std::vector<double>amp)
{
    int width, height, cx, cy;
    int num_v = amp.size();
    width = height = amp.size() * 2 + 10;
    cx = width / 2;
    cy = height / 2;
    cv::Mat img;
    img.create(width, height, 1);

    double amp_min = INF;
    double amp_max = -INF;

    for (std::vector<double>::iterator iter = amp.begin(); iter < amp.end(); iter++)
    {
        if (*iter < amp_min)
            amp_min = *iter;
        if (amp_min < *iter)
            amp_max = *iter;
    }
    double differ_dis = amp_max - amp_min;

    for (int i = 0; i < num_v; i++)
    {
        circle(img, cv::Point(cx, cy), i, (int)((amp[i]-amp_min)*255 / differ_dis));
    }

    if (img.rows > 600)
    {
        cv::Mat dst;
        double fx = (double)600 / img.cols;
        resize(img, dst,cv::Size(), fx, fx);
        imshow("PSF Plane", dst);
    }
    else
    {
        imshow("PSF Plane", img);
    }

    cv::waitKey(-1);

    return 0;
}

int simple_show_vec(const std::vector<std::vector<double> >&plane)
{
    if (plane.size() < 1) {
        std::cout << "Plane size not legal" << std::endl;
        return -1;
    }
    double max_pixel = plane[0][0];
    unsigned long height, width;
    height = plane.size();
    width = plane[0].size();
    cv::Mat img = cv::Mat::zeros(height, width, CV_8U);
    uchar *p;
    for (int i = 0; i < height; ++i) {
        p = img.ptr<uchar>(i);
        for (int j = 0; j < height; ++j) {
            p[j] = plane[i][j];
        }
    }
    cv::normalize(img, img, 255, 0);
    cv::resize(img, img, cv::Size(512, 512));
    cv::imshow("Direct Matrix", img);
    cv::waitKey(-1);
    return 0;
}

int show2DVec(const std::vector<std::vector<double> >&plane)
{
    if (plane.size() < 1)
    {
        std::cout << "Plane size not legal" << std::endl;
        return -1;
    }
    int width, height;
    height = plane.size();
    width = height;
    std::cout << "I got the Width: " << width << "\tHeight: " << height << std::endl;
    cv::Mat img = cv::Mat::zeros(height*2, width*2, CV_8U);
    uchar *p1, *p2;
    double max_pixel = plane[0][0];
    for (int i = 0; i < height; i++)
    {
        p1 = img.ptr<uchar>(height + i);
        p2 = img.ptr<uchar>(height - i);
        for (int j = 0; j < width; j++)
        {
            p1[width+j] = plane[i][j] * 255 / max_pixel;
            p1[width-j] = plane[i][j] * 255 / max_pixel;
            p2[width+j] = plane[i][j] * 255 / max_pixel;
            p2[width-j] = plane[i][j] * 255 / max_pixel;
        }
    }
    cv::imshow("Complete PSF", img);
    cv::imwrite("images/_Complete PSF.png", img);
    cv::waitKey(-1);
    return 0;
}

void vec2mat(std::vector<std::vector<double> > M2D, const cv::OutputArray out_im) {
    double *ptr;
    out_im.create(M2D.size(), M2D.size(), CV_64F);
    for (int i = 0; i < M2D.size(); ++i) {
        ptr = out_im.getMat().ptr<double>(i);
        for (int j = 0; j < M2D[0].size(); ++j) {
            ptr[j] = M2D[i][j];
        }
    }
}

int getTIFFinfo(const char* filename)
{
    TIFF* tif = TIFFOpen(filename, "r");

    if (tif) {
        int width, height, blacklevel;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
        TIFFGetField(tif, TIFFTAG_BLACKLEVEL, &blacklevel);

        std::cout << "TIFF_IMAGE info: "         << std::endl
                  << "WIDTH: "      << width     << "\t"
                  << "HEIGHT: "     << height    << std::endl
                  << "BLACK_LEVEL: "<< blacklevel<< std::endl;
        uint16 dircount = 0;
        do {
            dircount ++;
        } while (TIFFReadDirectory(tif));
        std::cout << "dircount: " << dircount << std::endl;
        TIFFClose(tif);
    }
}

int TIFFframenumber(const char* filename)
{
    TIFF* tif = TIFFOpen(filename, "r");
    int num;
    if (tif) {
        num = TIFFNumberOfDirectories(tif);
        TIFFClose(tif);
        return num;
    }
    return 0;
}

int getTIFF(const char* filename, cv::OutputArray out_img, uint16 frame_seq)
{
    // Open figure
    TIFF* tif = TIFFOpen(filename, "r");

    if (tif) {
        tdata_t buf;

        int width, height;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

        out_img.create(height, width, CV_64F);

        uint16 bitspersample = 1;
        uint16 samplesperpixel = 1;

        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitspersample);

        TIFFSetDirectory(tif, frame_seq);
        buf = _TIFFmalloc(TIFFScanlineSize(tif));
        double* ptr;
        for (int row = 0; row < height; ++row) {
            TIFFReadScanline(tif, buf, row);
            ptr = out_img.getMat().ptr<double>(row);
            for (int i = 0; i < width; ++i) {
                ptr[i] = *((uint16*)buf+i);
            }
        }
        _TIFFfree(buf);
    }
    TIFFClose(tif);
    return 0;
}

int reconstructZY(const char* filename, cv::OutputArray out_img)
{
    TIFF *tif = TIFFOpen(filename, "r");
    int num;
    if (tif) {
        tdata_t buf;
        int width, height;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

        out_img.create(width, height, CV_64F);
        buf = _TIFFmalloc(TIFFScanlineSize(tif));
        double *ptr;

        num = TIFFNumberOfDirectories(tif);
        // TODO(peo): Add highgui slide bar
        for (int i = 0; i < 0; ++i) {
            for (int row = 0; row < height; ++row) {

            }
        }
    }
}

static void on_trackbar(int , void*) {
    cv::Mat image;
    getTIFF(filename_TIFF, image, alpha_slider);
    cv::normalize(image, image, 255, 0);
    cv::imshow("NeuronXY", image);
}

int cv_slideWin(const char* filename)
{
    filename_TIFF = (char *) filename;
    int alpha_slider_max = TIFFframenumber(filename_TIFF);
    alpha_slider = 20;

    cv::namedWindow("NeuronXY");
    char TrackbarName[50];
    sprintf(TrackbarName, "Slice x %d", alpha_slider_max);
    cv::createTrackbar(TrackbarName, "NeuronXY", &alpha_slider, alpha_slider_max, on_trackbar);

    on_trackbar(alpha_slider, 0);

    cv::waitKey(0);
    return 0;
}