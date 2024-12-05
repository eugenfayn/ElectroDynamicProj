#include <iostream>
#include <complex>
#include "geometry/geometry.h"
#include "quadrature/quadrature.h"
#include "cmath"
#include <chrono>
#include <lapack.h>
#include <random>
#include <unistd.h>
#include <filesystem>


using namespace std;

namespace fs = std::filesystem;

namespace {
    double WAVENUMBER;

    void initializeGeometryType(const fs::path& inputFile) {
        std::string filename = inputFile.filename().string();
        if (filename.find("angle") != std::string::npos) {
            WAVENUMBER = 41.88790205;
        } else {
            WAVENUMBER = 62.83185307;
        }
        std::cout << "WAVENUMBER: " << WAVENUMBER << std::endl;
    }
}


void calcPovTok(Vertex &j, const int N, Triangle* &triangles, std::complex<double>* ves){
    Triangle currentTriangle;
    Vertex C;
    Vertex Xm;
    double wK[4]{-9./16, 25./48, 25./48, 25./48};
    double ksi4[4]{1./3, 3./5, 1./5, 1./5};
    double eta4[4]{1./3, 1./5, 3./5, 1./5};
    Vertex e_minus;
    Vertex e_plus;
    for (int i=0; i<N; i++){
        e_minus = Vertex(0,0,0,0);
        e_plus = Vertex(0,0,0,0);
        currentTriangle = triangles[2 * i];
        C = currentTriangle.getC();
        for (int m=0; m<4;m++){
            // p = 1, tri-
            Xm = currentTriangle.getA() * ksi4[m]+ currentTriangle.getB() *  eta4[m] +  currentTriangle.getC() * (1 - ksi4[m] - eta4[m]);
            e_minus = e_minus + (Xm - C)/currentTriangle.calc_S();
        }
        j = j + e_minus * ves[i];
        currentTriangle = triangles[2 * i + 1];
        C = currentTriangle.getC();
        for (int m=0; m<4;m++){
            // p = 2, tri+
            Xm = currentTriangle.getA() * ksi4[m]+ currentTriangle.getB() *  eta4[m] +  currentTriangle.getC() * (1 - ksi4[m] - eta4[m]);
            e_plus = e_plus + (C - Xm)/currentTriangle.calc_S();
        }
        j = j + e_plus * ves[i];
    }
}

// Generate random float number
double randomFloat() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(0.0, 1.0);
    return dis(gen);
}

// Добавить волновое число правильное
std::complex<double> calcF_Aij(Vertex Xk, Vertex Ym, Vertex Cx, Vertex Cy, double wavenumber=WAVENUMBER){
    std::complex<double> RES(0.0, 0.0);
    double R = Xk.distance(Ym);

    if (R < 1e-7) {
        return 0;
    }

    std::complex<double> exp_(cos(wavenumber*R), sin(wavenumber*R));

    Vertex Point1 = Cx - Xk; // Xk в квадратуре
    Vertex Point2 = Cy - Ym; // Ym в квадратуре

    double scalar = Point1.scalar_product(Point2);
    std::complex<double> POW(0.0, wavenumber * R);
    RES += std::complex<double>(0.0,0.0) + (-4.0 + wavenumber * wavenumber * scalar) * exp_/R;
    return RES;
}

// Достать список весов отдельно вынести
std::complex<double> FindI_Aij(Triangle sigmaX,Triangle sigmaY,Vertex Cx, Vertex Cy){
    double wK[4]{-9.0/16.0, 25.0/48.0, 25.0/48.0, 25.0/48.0}; // Четырёхточечный Гаусс. Веса.
    double wM[3]{1.0/3.0, 1.0/3.0, 1.0/3.0}; // Трёхточечный Гаусс. Веса.
    double ksi4[4]{1.0/3.0, 3.0/5.0, 1.0/5.0, 1.0/5.0};
    double eta4[4]{1.0/3.0, 1.0/5.0, 3.0/5.0, 1.0/5.0};
    double ksi3[3]{1.0/6.0, 2.0/3.0, 1.0/6.0};
    double eta3[3]{1.0/6.0, 1.0/6.0, 2.0/3.0};
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
void BuildMatrix(std::complex<double>** &A, const int N, const Triangle* &triangles){
    // Идём по I,J. 
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
            // p=1,q=2
            coef = - 1 / (4 * M_PI * triangles[2 * i].calc_S() * triangles[2 * j + 1].calc_S());
            sum_ += coef * FindI_Aij(triangles[2 * i],triangles[2 * j +1], triangles[2 * i].getC() , triangles[2 * j + 1].getC());
            // p=2,q=1
            coef = - 1 / (4 * M_PI * triangles[2 * i + 1].calc_S() * triangles[2 * j].calc_S());
            sum_ += coef * FindI_Aij(triangles[2 * i + 1],triangles[2 * j], triangles[2 * i + 1].getC() , triangles[2 * j].getC());
            // p=2,q=2
            coef = 1 / (4 * M_PI * triangles[2 * i + 1].calc_S() * triangles[2 * j + 1].calc_S());
            sum_ += coef * FindI_Aij(triangles[2 * i + 1],triangles[2 * j + 1], triangles[2 * i + 1].getC() , triangles[2 * j + 1].getC());
            // Записываем результат в A_ij
            A[i][j] = sum_;
        }
    }
}

// f - COMPLEX вектор правой части матричного уравнения для RWG
// N - размер вектора
// triangles - массив треугольников
// polarization - вектор поляризации, пока заговнокодим его координаты в Vertex
// tension- вектор напряжения, также
// wavenumber - волновое число
void BuildRightPart(std::complex<double>* &f, const int N, const Triangle* &triangles, const Vertex polarization, const Vertex tension, double wavenumber=WAVENUMBER){
    double wK[4]{-9./16, 25./48, 25./48, 25./48};
    // double wM[3]{1./3, 1./3, 1./3};
    double ksi4[4]{1./3, 3./5, 1./5, 1./5};
    double eta4[4]{1./3, 1./5, 3./5, 1./5};
    // double ksi3[3]{1./6, 2./3, 1./6};
    // double eta3[3]{1./6, 1./6, 2./3};
    for (int i=0; i < N; i++){
        std::complex<double> sum_(0.0,0.0);
        for (int p=0; p <2 ; p++){
            for (int k=0; k < 4 ; k++){
                int sign = (1)?(-1):(p%2);
                Triangle currentTriangle = triangles[2 * i + p%2];
                Vertex Cp = currentTriangle.getC();
                Vertex Xk = currentTriangle.getA() * ksi4[k]+ currentTriangle.getB() *  eta4[k] +  currentTriangle.getC() * (1 - ksi4[k] - eta4[k]) ; 
                double Sp = currentTriangle.calc_S();
                double tensionProductX = Xk.scalar_product(tension);
                double POW = wavenumber*tensionProductX;
                std::complex<double> exp_(cos(POW), sin(POW));
                sum_ += sign * ((Cp*(-1) + Xk)/Sp).scalar_product(polarization) * exp_ * wK[k]; //* ?
            }
        }
        f[i] = sum_;
    }
}

void SolveSLE(std::complex<double>** &A, const int N, std::complex<double>* &b){
    std::complex<double> *A_1d = new std::complex<double>[N * N];
    std::complex<double> *A_copy = new std::complex<double>[N * N];
    std::complex<double> *b_copy = (std::complex<double>*)malloc(N * sizeof(std::complex<double>));

    // Copy elements from 2D to 1D
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            A_1d[i * N + j] = A[j][i];
            A_copy[i * N + j] = A[j][i];
        }
    }

    std::cout << "\n2D Matrix:\n";
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    std::cout << "\n1D Matrix:\n";
    for(int i = 0; i < 9; i++) {
        std::cout << A_1d[i] << " ";
        if((i + 1) % 5 == 0) std::cout << "\n"; // New line every 5 elements for readability
    }
    std::cout << "\n";

    for (int i = 0; i < N; i++){
        b_copy[i] = b[i];
    }

    std::cout << "\nb:\n";
    for (int i = 0; i < 9; i++){
        std::cout << b[i] << " ";
    }
    
    __complex__ double* A_casted = reinterpret_cast<__complex__ double*>(A_1d);
    __complex__ double* B_casted = reinterpret_cast<__complex__ double*>(b);
	
    int Nrhs = 1;
    int *Ipvt = (int*)malloc(N * sizeof(N));
    int info;
    LAPACK_zgesv(&N, &Nrhs, A_casted, &N, Ipvt, B_casted, &N, &info);
    if (info != 0) {
        std::cout << "LAPACK solver failed with error code: " << info << std::endl;
        // Handle error
    }

    std::cout << "\nb after LAPACK:\n";
    for (int i = 0; i < 9; i++){
        std::cout << b[i] << " ";
    }

    double error = 0.0;
    // Compare
    std::complex<double>* rhs = new std::complex<double>[N];
    for (int i = 0; i < N; ++i) {
            rhs[i] = std::complex<double>(0.0,0.0);
            for (int j = 0; j < N; ++j) {
                    rhs[i] += A_copy[j * N + i] * b[j];
            }
            error += abs(b_copy[i] - rhs[i]);
        }

    std::cout << "\nFinal rhs vector:\n";
    for (int i = 0; i < 9; i++) {
        std::cout << rhs[i] << " ";
    }
    std::cout << "\n";

    delete[] rhs;

    printf("Rhs - b error %.20f", error);
    printf("\n\n");

    double relative_error = 0.0;
    double norm_b = 0.0;
    for(int i = 0; i < N; i++) {
        std::complex<double> sum(0.0, 0.0);
        for(int j = 0; j < N; j++) {
            sum += A[i][j] * b[j];
        }
        relative_error += std::abs(sum - b_copy[i]);
        norm_b += std::abs(b_copy[i]);
    }
    std::cout << "Relative error (rhs-b error divided by abs sum of b): " << relative_error/norm_b << std::endl;
}

// tau - единичный вектор, координаты храним в Vertex
// g - вектр b из SolveSLE  
// N - размер вектора
// triangles - массив треугольников
double calcEPR(Vertex tau, Vertex Gm, const int N, const Triangle* &triangles,  double wavenumber=WAVENUMBER){
    double wK[4]{-9./16, 25./48, 25./48, 25./48};
    double ksi4[4]{1./3, 3./5, 1./5, 1./5};
    double eta4[4]{1./3, 1./5, 3./5, 1./5};
    double sum_ = 0;
    for (int i=0; i < 2*N; i++){
        for (int m=0; m < 4 ; m++){
                Vertex Gm(1.0, 1.0, 1.0, -1);//я хз как его считать 
                Triangle currentTriangle = triangles[i];
                Vertex Xm = currentTriangle.getA() * ksi4[m]+ currentTriangle.getB() *  eta4[m] +  currentTriangle.getC() * (1 - ksi4[m] - eta4[m]) ; 
                double Si = currentTriangle.calc_S();
                double tauProductX = Xm.scalar_product(tau);
                double POW = -wavenumber*tauProductX;
                std::complex<double> exp_(cos(POW), sin(POW));
                Vertex V =  (exp_ * wavenumber * wavenumber * (Gm - tau*Gm.scalar_product(tau)) * wK[m] * Si);
                sum_ += 4*M_PI *(V.x*V.x + V.y*V.y + V.z*V.z);
        }
    }
    return sum_;
}