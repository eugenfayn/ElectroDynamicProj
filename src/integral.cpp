#include <iostream>
#include <complex.h>
#include "geometry/geometry.h"
#include "quadrature/quadrature.h"

// Добавить волновое число правильное
std::complex<double> calcF(Vertex Xk, Vertex Ym, Vertex Cx, Vertex Cy, double wavenumber=1){
    std::complex<double> RES(0.0, 0.0);
    double R = Xk.distance(Ym);
    Vertex Point1 = Cx - Xk; // Xk в квадратуре
    Vertex Point2 = Cy - Ym; // Ym в квадратуре
    double scalar = Point1.scalar_product(Point2);
    std::complex<double> POW(0.0, wavenumber * R);
    RES += std::complex<double>(-4.0,0.0) + wavenumber * wavenumber * scalar * exp(POW)/R;
    return RES;
}

// Достать список весов отдельно вынести
std::complex<double> FindI(Triangle sigmaX,Triangle sigmaY,Vertex Cx, Vertex Cy){
    double wK[4]{-9/16, 25/48, 25/48, 25/48}; // Четырёхточечный Гаусс. Веса.
    double wM[3]{1/3, 1/3, 1/3}; // Трёхточечный Гаусс. Веса.
    std::complex<double> Xk;
    std::complex<double> Ym;
    std::complex<double> RES = sigmaX.calcSquare() * sigmaY.calcSquare();
    std::complex<double> SUM_(0.0,0.0);
    // Формирование интеграла
    for (int k=0;k<4;k++){
        for (int m=0;m<3;m++){
            SUM_ = SUM_ + wK[k] * wM[m] * calcF(Xk,Ym,Cx,Cy);
        }
    }
    RES = RES * SUM_;
    return RES;
};


// spisok Edges [0, 1, 2, 3]
// spisok Triangles [0-, 0+, 1-,1+, ]
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


