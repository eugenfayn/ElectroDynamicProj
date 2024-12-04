#include <iostream>
#include <complex.h>
#include "geometry/geometry.h"
#include "quadrature/quadrature.h"
#include "cmath"

// Добавить волновое число правильное
std::complex<double> calcF_Aij(Vertex Xk, Vertex Ym, Vertex Cx, Vertex Cy, double wavenumber=1){
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
std::complex<double> FindI_Aij(Triangle sigmaX,Triangle sigmaY,Vertex Cx, Vertex Cy){
    double wK[4]{-9/16, 25/48, 25/48, 25/48}; // Четырёхточечный Гаусс. Веса.
    double wM[3]{1/3, 1/3, 1/3}; // Трёхточечный Гаусс. Веса.
    double ksi4[4]{1/3, 3/5, 1/5, 1/5};
    double eta4[4]{1/3, 1/5, 3/5, 1/5};
    double ksi3[3]{1/6, 2/3, 1/6};
    double eta3[3]{1/6, 1/6, 2/3};
    std::complex<double> RES = sigmaX.calcSquare() * sigmaY.calcSquare();
    std::complex<double> SUM_(0.0,0.0);
    // Формирование интеграла
    for (int k=0;k<4;k++){
        Vertex Xk = sigmaX.a * ksi4[k]+ sigmaX.b *  eta4[k] +  sigmaX.c * (1 - ksi4[k] - eta4[k]) ; // Четырёхточечный Гаусс. Точка.
        for (int m=0;m<3;m++){
            Vertex Ym = sigmaY.a * ksi3[m]+ sigmaY.b *  eta3[m] +  sigmaY.c * (1 - ksi3[m] - eta3[m]); // Четырёхточечный Гаусс. Точка.
            SUM_ = SUM_ + wK[k] * wM[m] * calcF_Aij(Xk,Ym,Cx,Cy);
        }
    }
    RES = RES * SUM_;
    return RES;
};

// A - COMPLEX матрица, куда будут записаны коэффициенты A_ij
// NxN - размер матриц
// Чекнуть нумерацию рёбер с нуля или нет. В Edge начинается с нуля
// triangles - массив треугольников
void BuildMatrix(std::complex<double>** &A, int N, int b, Triangle* &trinagles){
    // Идём по I,J.
    // В случае если I,J равны 0    
    double coef = 0.0;
    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            std::complex<double> sum_(0.0,0.0);
            // p=1,q=1
            coef = 1 / (4 * M_PI );
            sum_ += 0;
            // p=1,q=2
            sum_ += 0;
            // p=2,q=1
            sum_ += 0;
            // p=2,q=2
            sum_ += 0;
            // Записываем результат в A_ij
            A[i][j] = sum_;
        }
    }
}


