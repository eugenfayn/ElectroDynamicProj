#include <iostream>
#include <complex.h>
#include "geometry/geometry.h"
#include "quadrature/quadrature.h"

std::complex<double> calcF(Vertex Xk, Vertex Ym, Vertex Cx, Vertex Cy, double wavenumber){
    std::complex<double> RES(0.0, 0.0);
    double R = Xk.abs(Ym);
    std::complex<double> POW(0)
    RES = -4 + wavenumber * wavenumber * (Cx - Xk) * (Cy -Ym) * exp(-)
    return RES;
}
// Достать список весов отдельно вынести
std::complex<double> FindI(Triangle sigmaX,Triangle sigmaY,Vertex Cx, Vertex Cy){
    double wK[4]{-9/16, 25/48, 25/48, 25/48}; // Четырёхточечный Гаусс. Веса.
    double wM[3]{1/3, 1/3, 1/3}; // Трёхточечный Гаусс. Веса.
    std::complex<double> RES = sigmaX.calcSquare() * sigmaY.calcSquare();
    std::complex<double> SUM_(0.0,0.0);
    // Формирование интеграла
    for (int k=0;k<4;k++){
        for (int m=0;m<3;m++){
            f

        }
    }
};


// spisok Edges [0,1,2,3]
// spisok Triangles [0-,0+,1-,1+]
//


// A - COMPLEX матрица, куда будут записаны коэффициенты A_ij
// NxN - размер матрицы
// Чекнуть нумерацию рёбер с нуля или нет. В Edge начинается с нуля
void BuildMatrix(std::complex<double>** &A, int N, int b){
    // Идём по I,J.
    // В случае если I,J равны 0    
    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            std::complex<double> sum_(0.0,0.0);
            // sigmaI- sigmaJ-
            sum_ += 0;
            // sigmaI- sigmaJ+
            sum_ += 0;
            // sigmaI+ sigmaJ-
            sum_ += 0;
            // sigmaI+ sigmaJ+
            sum_ += 0;
            // Записываем результат в A_ij
            A[i][j] = sum_;
        }
    }
}


