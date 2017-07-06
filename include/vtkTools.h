//
// Created by peo on 17-6-27.
//

#ifndef PSF_VTKTOOLS_H
#define PSF_VTKTOOLS_H

/* make some visualization tools for data showing
 * Author: hzn
 */
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <complex>
#include <cmath>

#include <vtk-7.1/vtkVersion.h>
#include <vtk-7.1/vtkRenderer.h>
#include <vtk-7.1/vtkRenderWindowInteractor.h>
#include <vtk-7.1/vtkSmartPointer.h>
#include <vtk-7.1/vtkChartXY.h>
#include <vtk-7.1/vtkTable.h>
#include <vtk-7.1/vtkPlot.h>
#include <vtk-7.1/vtkDoubleArray.h>
#include <vtk-7.1/vtkContextView.h>
#include <vtk-7.1/vtkContextScene.h>
#include <vtk-7.1/vtkPen.h>

/* 绘制二维图像的工具
 * 传入的都是双精度类型
 */
int vtk_2Dplot(std::vector<double> X, std::vector<double>Y);

#endif //PSF_VTKTOOLS_H
