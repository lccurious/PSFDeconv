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

#include <vtk-7.1/vtkVersion.h>
#include <vtk-7.1/vtkRenderer.h>
#include <vtk-7.1/vtkRenderWindow.h>
#include <vtk-7.1/vtkRenderWindowInteractor.h>
#include <vtk-7.1/vtkSmartPointer.h>
#include <vtk-7.1/vtkChartXY.h>
#include <vtk-7.1/vtkTable.h>
#include <vtk-7.1/vtkPlot.h>
#include <vtk-7.1/vtkChartLegend.h>
#include <vtk-7.1/vtkDoubleArray.h>
#include <vtk-7.1/vtkContextView.h>
#include <vtk-7.1/vtkContextScene.h>
#include <vtk-7.1/vtkPen.h>
#include <vtk-7.1/vtkPNGWriter.h>
#include <vtk-7.1/vtkWindowToImageFilter.h>

/* 绘制二维图像的工具
 * 传入的都是双精度类型
 */
int vtk_2Dplot(std::vector<double> X, std::vector<double>Y);
int vtk_2Dplot_com(std::vector<double> Y1, std::vector<double> Y2,
                   const char *strY1, const char *strY2);

#endif //PSF_VTKUTILS_H
