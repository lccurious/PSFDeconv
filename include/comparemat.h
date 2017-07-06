//
// Created by peo on 17-7-6.
//

#ifndef PSF_COMPAREMAT_H
#define PSF_COMPAREMAT_H
#include <iostream>
#include <string>
#include <vector>
#include <matio.h>

int read_mat(const char *filename);

int mat2vector(const char *filename, std::vector<std::vector<double> >& M2D);
int mat3vector(const char *filename,
               std::vector<std::vector<std::vector<double> > >& M3D);

#endif //PSF_COMPAREMAT_H
