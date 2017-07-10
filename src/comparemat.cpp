//
// Created by peo on 17-7-6.
//

#include "comparemat.h"

int read_mat(const char* filename)
{
    mat_t *mat = Mat_Open(filename, MAT_ACC_RDONLY);
    unsigned int indice = 0;
    if (mat) {
        std::cout << "mat open successfully" << std::endl;

        matvar_t *matVar = 0;
        // 验证常规级别的BW模型生成数据
        matVar = Mat_VarRead(mat, (char*)"psf");
        if (matVar) {
            unsigned xSize = matVar->nbytes/matVar->data_size;
            std::cout<< xSize << std::endl;
            std::cout<< xSize/matVar->dims[0] << std::endl;
            std::cout<< xSize/matVar->dims[1]/matVar->dims[0] << std::endl;
            const double *xData = static_cast<const double*>(matVar->data);
            if (matVar->rank == 3) {
                // loop1: stack depth
//                for (int i = 0; i < matVar->dims[2]; ++i) {
                    // loop2: psf width
                    for (int j = 0; j < matVar->dims[1]; ++j) {
                        // loop3: psf height
                        for (int k = 0; k < matVar->dims[0]; ++k) {
                            std::cout<<xData[indice+k]
                                     <<" ";
                        }
                        std::cout<<std::endl;
                        indice += matVar->dims[1];
                    }
                    std::cout<<"\n===================\n";
//                    indice += matVar->dims[0]*matVar->dims[1];
//                }
            }
        }
        Mat_Close(mat);
    }
    return 0;
}

int mat2vector(const char *filename, std::vector<std::vector<double> >& M2D,
               int stack_seq)
{
    mat_t *mat = Mat_Open(filename, MAT_ACC_RDONLY);
    unsigned int indice = 0;
    if (mat) {

        matvar_t *matVar = 0;
        // 验证常规级别的BW模型生成数据
        matVar = Mat_VarRead(mat, (char*)"psf_best");
        if (matVar) {
            std::cout << "mat open successfully" << std::endl;
            indice += matVar->dims[0]*matVar->dims[1]*stack_seq;
            const double *xData = static_cast<const double*>(matVar->data);
            if (matVar->rank == 3) {
                // loop1: psf rows
                M2D.resize(matVar->dims[1]);
                for (int j = 0; j < matVar->dims[1]; ++j) {
                    // loop2: psf columns
                    for (int k = 0; k < matVar->dims[0]; ++k) {
                        M2D[j].push_back(xData[indice+k]);
                    }
                    indice += matVar->dims[1];
                }
            }
        }
        std::cout<<"feed into vector finished" << std::endl;
        Mat_Close(mat);
    }
    return 0;
}

int mat3vector(const char *filename,
               std::vector<std::vector<std::vector<double> > >& M3D)
{
    mat_t *mat = Mat_Open(filename, MAT_ACC_RDONLY);
    unsigned int indice = 0;
    if (mat) {
        std::cout << "mat open successfully" << std::endl;

        matvar_t *matVar = 0;
        // 验证常规级别的BW模型生成数据
        matVar = Mat_VarRead(mat, (char*)"psf");
        if (matVar) {
            const double *xData = static_cast<const double*>(matVar->data);
            if (matVar->rank == 3) {
                M3D.resize(matVar->dims[2]);
                for (int i=0; i < matVar->dims[2]; ++i) {
                    // loop1: psf rows
                    M3D[i].resize(matVar->dims[1]);
                    for (int j = 0; j < matVar->dims[1]; ++j) {
                        // loop2: psf columns
                        for (int k = 0; k < matVar->dims[0]; ++k) {
                            M3D[i][j].push_back(xData[indice+k]);
                        }
                        indice += matVar->dims[1];
                    }
                }
            }
        }
        std::cout<<"feed into vector finished" << std::endl;
        Mat_Close(mat);
    }
    return 0;
}

int vtk_compare_linear(const char *filename,
                       std::vector<double> Y1,
                       std::vector<double> Y2)
{
    assert(Y1.size() == Y2.size());
    for (int i = 0; i < Y1.size(); ++i) {
        std::cout << Y1[i] << "\t" << Y2[i] << std::endl;
    }
    return 0;
}