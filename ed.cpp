#include <iostream>
#include <complex>

using namespace std;
// Чтение геометрии
//  интегрирования произвольной функции
// Calculating Aij, Bi
// Решение системы
// Постобработка

// func for calculating A_ij

extern "C"
{
	extern void zgesv_(int *N, int *Nrhs, complex<double> *A, int *ldA ,     int *Ipvt, complex<double> *B, int *ldB, int *info);
}

int main(){
	cout << "Geometry" << endl;
	return 0;
}
