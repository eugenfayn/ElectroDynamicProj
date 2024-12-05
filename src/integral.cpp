#include <iostream>
#include <complex.h>
#include "geometry/geometry.h"
#include "quadrature/quadrature.h"
#include "cmath"
// добавить eps=10-6
// Добавить волновое число правильное
std::complex<double> calcF_Aij(Vertex Xk, Vertex Ym, Vertex Cx, Vertex Cy, double wavenumber=1){
    std::complex<double> RES(0.0, 0.0);
    double R = Xk.distance(Ym);
    if (R<1e-6){
        return 0;
    }
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
    double ksi4[4]{1/3, 3/5, 1/5, 1/5}; // ksi четырёхточчный Гаусс.
    double eta4[4]{1/3, 1/5, 3/5, 1/5}; // eta четырёхточечны Гаусс.
    double ksi3[3]{1/6, 2/3, 1/6}; // ksi трёхточечный Гаусс.
    double eta3[3]{1/6, 1/6, 2/3}; // eta трёхточечны Гаусс.
    std::complex<double> RES = sigmaX.calc_S() * sigmaY.calc_S();
    std::complex<double> SUM_(0.0,0.0);
    // Формирование интеграла
    for (int k=0;k<4;k++){
        Vertex Xk = sigmaX.getA() * ksi4[k] + sigmaX.getB() *  eta4[k] +  sigmaX.getC() * (1 - ksi4[k] - eta4[k]) ; // Четырёхточечный Гаусс. Точка.
        for (int m=0;m<3;m++){
            Vertex Ym = sigmaY.getA() * ksi3[m]+ sigmaY.getB() *  eta3[m] +  sigmaY.getC() * (1 - ksi3[m] - eta3[m]); // Четырёхточечный Гаусс. Точка.
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
void BuildMatrix(std::complex<double>** &A, int N, int b, Triangle* &triangles){
    // Идём по I,J.
    // В случае если I,J равны 0    
    double coef = 0.0;
    for (int i=0; i < N; i++){
        for (int j=0; j < N; j++){
            if (i==j){
                A[i][j] = std::complex<double>(0.0,0.0);
            }
            else{
            std::complex<double> sum_(0.0,0.0);
            // p=1,q=1
            coef = 1 / (4 * M_PI * triangles[2 * i].calc_S() * triangles[2 * j].calc_S());
            sum_ += coef * FindI_Aij(triangles[2 * i],triangles[2 * j], triangles[2 * i].getC() , triangles[2 * j].getC());
            std::cout << "SUM 1: " << sum_ << std::endl;
            std::cout << "coef 1 " << coef << std::endl;
            // p=1,q=2
            coef = - 1 / (4 * M_PI * triangles[2 * i].calc_S() * triangles[2 * j + 1].calc_S());
            sum_ += coef * FindI_Aij(triangles[2 * i],triangles[2 * j +1], triangles[2 * i].getC() , triangles[2 * j + 1].getC());
            std::cout << "SUM 2: " << sum_ << std::endl;
            std::cout << "coef 2 " << coef << std::endl;
            // p=2,q=1
            coef = - 1 / (4 * M_PI * triangles[2 * i + 1].calc_S() * triangles[2 * j].calc_S());
            sum_ += coef * FindI_Aij(triangles[2 * i + 1],triangles[2 * j], triangles[2 * i + 1].getC() , triangles[2 * j].getC());
            std::cout << "SUM 3: " << sum_ << std::endl;
            std::cout << "coef 3 " << coef << std::endl;
            // p=2,q=2
            coef = 1 / (4 * M_PI * triangles[2 * i + 1].calc_S() * triangles[2 * j + 1].calc_S());
            sum_ += coef * FindI_Aij(triangles[2 * i + 1],triangles[2 * j + 1], triangles[2 * i + 1].getC() , triangles[2 * j + 1].getC());
            std::cout << "SUM 4: " << sum_ << std::endl;
            std::cout << "coef 4 " << coef << std::endl;
            // Записываем результат в A_ij
            A[i][j] = sum_;
            }
        }
    }
}


// tau - единичный вектор, координаты храним в Vertex
// g - вектр b из SolveSLE  
// N - размер вектора
// triangles - массив треугольников
double calcEPR(Vertex tau, std::complex<double>* &g, const int N, const Triangle* &triangles,  double wavenumber=1){
    double wK[4]{-9./16, 25./48, 25./48, 25./48};
    double ksi4[4]{1./3, 3./5, 1./5, 1./5};
    double eta4[4]{1./3, 1./5, 3./5, 1./5};
    double sum_ = 0;
    double* G[N];
    for (int i=0; i < 2*N; i++){
        for (int m=0; m < 4 ; m++){
                Vertex Gm(1.0, 1.0, 1.0, -1);//я хз как его считать
                Triangle currentTriangle = triangles[i];
                Vertex Xm = currentTriangle.getA() * ksi4[m]+ currentTriangle.getB() *  eta4[m] +  currentTriangle.getC() * (1 - ksi4[m] - eta4[m]) ; 
                double Si = currentTriangle.calc_S();
                double tauProductX = Xm.scalar_product(tau);
                double POW = -wavenumber*tauProductX;
                std::complex<double> exp_(cos(POW), sin(POW));
                Vertex V =  (exp_ * wavenumber * wavenumber * (G[m] - tau*Gm.scalar_product(tau)) * wK[m] * Si);
                sum_ += 4*M_PI *(V.x*V.x + V.y*V.y + V.z*V.z);
        }
    }
    return sum_;
}