//
// Created by peo on 17-7-10.
//

#include "vtkutils.h"

int vtk_2Dplot(std::vector<double> X, std::vector<double> Y)
{
    if (X.size() != Y.size())
    {
        std::cout << "X Y length not equal!" << std::endl;
        return -1;
    }

    vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkDoubleArray> arrX =
            vtkSmartPointer<vtkDoubleArray>::New();

    // TODO(peo):Modify into matplotlib style
    arrX->SetName("X Axis");
    table->AddColumn(arrX);

    vtkSmartPointer<vtkDoubleArray>arrC =
            vtkSmartPointer<vtkDoubleArray>::New();
    arrC->SetName("Function Value");
    table->AddColumn(arrC);

    int num_point = X.size();
    table->SetNumberOfRows(num_point);
    for (int i = 0; i < num_point; i++)
    {
        table->SetValue(i, 0, X[i]);
        table->SetValue(i, 1, Y[i]);
    }

    // Setup the View
    vtkSmartPointer<vtkContextView>view =
            vtkSmartPointer<vtkContextView>::New();
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    // Add line plots, setting the colors etc;
    vtkSmartPointer<vtkChartXY> chart =
            vtkSmartPointer<vtkChartXY>::New();
    view->GetScene()->AddItem(chart);
    vtkPlot *line = chart->AddPlot(vtkChart::LINE);

    line->SetInputData(table, 0, 1);
    line->SetColor(0, 255, 0, 255);
    line->SetWidth(2.0);
    line = chart->AddPlot(vtkChart::LINE);

    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();

    return EXIT_SUCCESS;
}

int vtk_2Dplot_com(std::vector<double> Y1, std::vector<double> Y2)
{
    if (Y2.size() != Y1.size())
    {
        std::cout << "X Y length not equal!" << std::endl;
        return -1;
    }

    vtkSmartPointer<vtkTable> table =
            vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkDoubleArray> arrX =
            vtkSmartPointer<vtkDoubleArray>::New();

    // TODO(peo):Modify into matplotlib style
    arrX->SetName("X Axis");
    table->AddColumn(arrX);

    vtkSmartPointer<vtkDoubleArray> arrY1 =
            vtkSmartPointer<vtkDoubleArray>::New();
    arrY1->SetName("Y1 Value");
    table->AddColumn(arrY1);

    vtkSmartPointer<vtkDoubleArray> arrY2 =
            vtkSmartPointer<vtkDoubleArray>::New();
    arrY2->SetName("Y2 Value");
    table->AddColumn(arrY2);

    int num_point = Y1.size();
    table->SetNumberOfRows(num_point);
    for (int i = 0; i < num_point; i++)
    {
        table->SetValue(i, 0, (double)i);
        table->SetValue(i, 1, Y1[i]);
        table->SetValue(i, 2, Y2[i]);
    }

    // Setup the View
    vtkSmartPointer<vtkContextView>view =
            vtkSmartPointer<vtkContextView>::New();
    view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

    // Add line plots, setting the colors etc;
    vtkSmartPointer<vtkChartXY> chart =
            vtkSmartPointer<vtkChartXY>::New();
    view->GetScene()->AddItem(chart);
    vtkPlot *line = chart->AddPlot(vtkChart::LINE);

    line->SetInputData(table, 0, 1);
    line->SetColor(255, 0, 0, 255);
    line->SetWidth(2.0);

    line = chart->AddPlot(vtkChart::LINE);
    line->SetInputData(table, 0, 2);
    line->SetColor(0, 0, 255, 255);
    line->SetWidth(2.0);

    view->GetInteractor()->Initialize();
    view->GetInteractor()->Start();

    return EXIT_SUCCESS;
}
