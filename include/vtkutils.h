//
// Created by peo on 17-7-10.
//

#ifndef PSF_VTKUTILS_H
#define PSF_VTKUTILS_H


#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <complex>
#include <cmath>

#include <vtk-8.0/vtkVersion.h>
#include <vtk-8.0/vtkRenderer.h>
#include <vtk-8.0/vtkRenderWindow.h>
#include <vtk-8.0/vtkRenderWindowInteractor.h>
#include <vtk-8.0/vtkSmartPointer.h>
#include <vtk-8.0/vtkChartXY.h>
#include <vtk-8.0/vtkTable.h>
#include <vtk-8.0/vtkPlot.h>
#include <vtk-8.0/vtkChartLegend.h>
#include <vtk-8.0/vtkDoubleArray.h>
#include <vtk-8.0/vtkContextView.h>
#include <vtk-8.0/vtkContextScene.h>
#include <vtk-8.0/vtkPen.h>
#include <vtk-8.0/vtkPNGWriter.h>
#include <vtk-8.0/vtkWindowToImageFilter.h>

/* 绘制二维图像的工具
 * 传入的都是双精度类型
 */
int vtk_2Dplot(std::vector<double> X, std::vector<double>Y);
int vtk_2Dplot_com(std::vector<double> Y1, std::vector<double> Y2,
                   const char *strY1, const char *strY2);
int vtk1Dplot(const std::vector<double> Y);

#endif //PSF_VTKUTILS_H
